//! Short-Weierstrass arithmetic for Bandersnatch and the paper's §5.3
//! exponentiation gadget (chord rule, 3 constraints per bit).
//!
//! Twisted Edwards  a·x² + y² = 1 + d·x²y²
//!   → Montgomery   B_m·v² = u³ + A_m·u² + u        (u = (1+y)/(1-y), v = u/x)
//!   → Weierstrass  Y² = X³ + A·X + B               (X = u/B_m + A_m/(3B_m), Y = v/B_m)
//!
//! Written with AI assistance (Claude), following Golden §5.3.

use ark_ec::twisted_edwards::TECurveConfig;
use ark_ec::{AffineRepr, CurveGroup};
use ark_ed_on_bls12_381_bandersnatch::{
    BandersnatchConfig, EdwardsAffine, EdwardsProjective, Fq, Fr as Scalar,
};
use ark_ff::{Field, One, PrimeField, Zero};
use ark_relations::{
    lc,
    r1cs::{ConstraintSystemRef, SynthesisError, Variable},
};

// ───────────────────────── W1: conversion and native arithmetic ─────────────────────────

/// Weierstrass constants (A, B) and the Montgomery pair (A_m, B_m) for Edwards (a, d).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SwParams<F> {
    pub a_sw: F,
    pub b_sw: F,
    pub a_m: F,
    pub b_m: F,
}

pub fn sw_params<F: PrimeField>(a: F, d: F) -> SwParams<F> {
    let a_minus_d_inv = (a - d).inverse().expect("a != d");
    let a_m = F::from(2u64) * (a + d) * a_minus_d_inv;
    let b_m = F::from(4u64) * a_minus_d_inv;
    let three = F::from(3u64);
    let a_sw = (three - a_m.square()) * (three * b_m.square()).inverse().unwrap();
    let b_sw = (F::from(2u64) * a_m.square() * a_m - F::from(9u64) * a_m)
        * (F::from(27u64) * b_m.square() * b_m).inverse().unwrap();
    SwParams { a_sw, b_sw, a_m, b_m }
}

/// The Bandersnatch parameters, computed once.
pub fn bandersnatch_sw() -> SwParams<Fq> {
    sw_params(BandersnatchConfig::COEFF_A, BandersnatchConfig::COEFF_D)
}

/// Edwards (x, y) → Weierstrass (X, Y). Not defined for the identity (0, 1)
/// (it maps to the point at infinity) or for (0, -1).
pub fn te_to_sw<F: PrimeField>(p: &SwParams<F>, x: F, y: F) -> (F, F) {
    let u = (F::one() + y) * (F::one() - y).inverse().expect("identity has no affine SW form");
    let v = u * x.inverse().expect("x = 0 has no affine SW form");
    let b_m_inv = p.b_m.inverse().unwrap();
    let xx = u * b_m_inv + p.a_m * (F::from(3u64) * p.b_m).inverse().unwrap();
    let yy = v * b_m_inv;
    (xx, yy)
}

/// Convenience: an Edwards affine point → Weierstrass coordinates.
pub fn point_to_sw(p: &EdwardsAffine) -> (Fq, Fq) {
    te_to_sw(&bandersnatch_sw(), p.x, p.y)
}

/// Y² == X³ + A·X + B
pub fn on_sw_curve<F: PrimeField>(p: &SwParams<F>, x: F, y: F) -> bool {
    y.square() == x.square() * x + p.a_sw * x + p.b_sw
}

/// Chord rule, P + Q with x1 != x2.
pub fn sw_add_native<F: PrimeField>(x1: F, y1: F, x2: F, y2: F) -> (F, F) {
    let s = (y2 - y1) * (x2 - x1).inverse().expect("chord rule needs x1 != x2");
    let x3 = s.square() - x1 - x2;
    let y3 = s * (x1 - x3) - y1;
    (x3, y3)
}

/// Tangent rule, 2P (y != 0).
pub fn sw_double_native<F: PrimeField>(p: &SwParams<F>, x: F, y: F) -> (F, F) {
    let s = (F::from(3u64) * x.square() + p.a_sw) * (F::from(2u64) * y).inverse().unwrap();
    let x3 = s.square() - x - x;
    let y3 = s * (x - x3) - y;
    (x3, y3)
}

// ───────────────────────── W2: precomputed tables and the gadget ─────────────────────────

/// For each bit j: (C_j, D_j) in Weierstrass coordinates, where
///   P_j = 2^j · base,   C_j = c_j · h_c,   D_j = P_j + C_j,
/// and Δ_j = D_j if bit j is 1, C_j if it is 0.
///
/// Correction exponents (paper §5.3, with the typo fixed):
///   c_0 = 1,  c_j = 2 (0 < j < λ),  c_λ = s - 2λ + 1  ≡ -(2λ - 1),
/// so Σ c_j ≡ 0 (mod s) and the corrections cancel: L_λ = k · base.
/// The corrections use an independent generator h_c (not the base) so the
/// partial sums can't be steered onto equal x-coordinates.
pub fn sw_fixed_table(
    base: &EdwardsAffine,
    h_c: &EdwardsAffine,
    nbits: usize,
) -> Vec<((Fq, Fq), (Fq, Fq))> {
    assert!(nbits >= 2);
    let lambda = nbits - 1;
    let swp = bandersnatch_sw();
    let base_p: EdwardsProjective = (*base).into();
    let hc_p: EdwardsProjective = (*h_c).into();

    let mut table = Vec::with_capacity(nbits);
    let mut p_j = base_p; // 2^j · base
    for j in 0..nbits {
        let c_j: EdwardsProjective = if j == 0 {
            hc_p
        } else if j < lambda {
            hc_p + hc_p
        } else {
            // -(2λ - 1) · h_c, the same as (s - 2λ + 1) · h_c
            -(hc_p * Scalar::from((2 * lambda - 1) as u64))
        };
        let d_j = p_j + c_j;
        let c_a = c_j.into_affine();
        let d_a = d_j.into_affine();
        assert!(!c_a.is_zero() && !d_a.is_zero(), "table point at infinity (negligible)");
        table.push((te_to_sw(&swp, c_a.x, c_a.y), te_to_sw(&swp, d_a.x, d_a.y)));
        p_j = p_j + p_j;
    }
    table
}

/// k · base for a PUBLIC base, via precomputed (C_j, D_j) and the chord rule.
/// `bits` LSB first. Cost: 2 + 3·(nbits - 1) constraints.
/// Returns Weierstrass coordinates of k · base.
pub fn scalar_mul_sw_fixed(
    cs: ConstraintSystemRef<Fq>,
    bits: &[Variable],
    bit_vals: &[Fq],
    table: &[((Fq, Fq), (Fq, Fq))],
) -> Result<(Variable, Variable, Fq, Fq), SynthesisError> {
    assert_eq!(bits.len(), table.len());
    assert_eq!(bit_vals.len(), table.len());

    // Δ_j as a linear combination of the bit:  C + k·(D - C)
    let delta_lc = |j: usize| {
        let ((xc, yc), (xd, yd)) = table[j];
        let x = lc!() + (xc, Variable::One) + (xd - xc, bits[j]);
        let y = lc!() + (yc, Variable::One) + (yd - yc, bits[j]);
        (x, y)
    };
    let delta_val = |j: usize| {
        let ((xc, yc), (xd, yd)) = table[j];
        if bit_vals[j].is_one() { (xd, yd) } else { (xc, yc) }
    };

    // Base case: L_0 = Δ_0  (2 constraints)
    let (mut lx, mut ly) = delta_val(0);
    let mut lx_var = cs.new_witness_variable(|| Ok(lx))?;
    let mut ly_var = cs.new_witness_variable(|| Ok(ly))?;
    let (dx0, dy0) = delta_lc(0);
    cs.enforce_constraint(lc!() + lx_var, lc!() + Variable::One, dx0)?;
    cs.enforce_constraint(lc!() + ly_var, lc!() + Variable::One, dy0)?;

    // L_j = L_{j-1} + Δ_j  (3 constraints each, chord rule)
    for j in 1..table.len() {
        let (dxv, dyv) = delta_val(j);
        let denom = lx - dxv;
        if denom.is_zero() {
            // equal x-coordinates: excluded by the corrections except with negligible probability
            return Err(SynthesisError::Unsatisfiable);
        }
        let s = (ly - dyv) * denom.inverse().unwrap();
        let nx = s.square() - lx - dxv;
        let ny = s * (lx - nx) - ly;

        let s_var = cs.new_witness_variable(|| Ok(s))?;
        let nx_var = cs.new_witness_variable(|| Ok(nx))?;
        let ny_var = cs.new_witness_variable(|| Ok(ny))?;
        let (dx, dy) = delta_lc(j);

        // s · (x_L − x_Δ) = y_L − y_Δ
        cs.enforce_constraint(lc!() + s_var, lc!() + lx_var - &dx, lc!() + ly_var - &dy)?;
        // s · s = x_L + x_L' + x_Δ
        cs.enforce_constraint(lc!() + s_var, lc!() + s_var, lc!() + lx_var + nx_var + &dx)?;
        // s · (x_L − x_L') = y_L + y_L'
        cs.enforce_constraint(lc!() + s_var, lc!() + lx_var - nx_var, lc!() + ly_var + ny_var)?;

        lx = nx;
        ly = ny;
        lx_var = nx_var;
        ly_var = ny_var;
    }
    Ok((lx_var, ly_var, lx, ly))
}
