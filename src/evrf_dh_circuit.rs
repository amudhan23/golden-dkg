use ark_ff::{BigInteger, PrimeField};
use ark_relations::{
    lc,
    r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError, Variable},
};

use crate::{
    circuit_evrf::bit_decompose,
    scalar_mul::{scalar_mul, scalar_mul_fixed},
};

pub fn evrf_dh<F: PrimeField>(
    cs: ConstraintSystemRef<F>,
    g_in_x: F,
    g_in_y: F,
    pk1_x: F,
    pk1_y: F,
    pk2_x: F,
    pk2_y: F,
    h1_x: F,
    h1_y: F,
    h2_x: F,
    h2_y: F,
    beta: F,
    witness: F,
    a: F,
    d: F,
) -> Result<(Variable, F), SynthesisError> {
    let sk_var = cs.new_witness_variable(|| Ok(witness))?;

    let (sk_bits_var, sk_bits) = bit_decompose(cs.clone(), witness, sk_var)?;

    let (pk1_cons_x_var, pk1_cons_y_var, _, _) =
        scalar_mul_fixed(cs.clone(), &sk_bits_var, &sk_bits, g_in_x, g_in_y, a, d)?;

    let pk1_x_var = cs.new_input_variable(|| Ok(pk1_x))?;
    let pk1_y_var = cs.new_input_variable(|| Ok(pk1_y))?;

    cs.enforce_constraint(
        lc!() + pk1_x_var,
        lc!() + Variable::One,
        lc!() + pk1_cons_x_var,
    )?;

    cs.enforce_constraint(
        lc!() + pk1_y_var,
        lc!() + Variable::One,
        lc!() + pk1_cons_y_var,
    )?;

    let (s_x_var, _, s_x, _) =
        scalar_mul_fixed(cs.clone(), &sk_bits_var, &sk_bits, pk2_x, pk2_y, a, d)?;

    let (k_bits_var, k_bits) = bit_decompose(cs.clone(), s_x, s_x_var)?;

    let (t1_x_var, _, t1_x, _) =
        scalar_mul_fixed(cs.clone(), &k_bits_var, &k_bits, h1_x, h1_y, a, d)?;

    let (t2_x_var, _, t2_x, _) =
        scalar_mul_fixed(cs.clone(), &k_bits_var, &k_bits, h2_x, h2_y, a, d)?;

    let r = beta * t1_x + t2_x;

    let r_var = cs.new_witness_variable(|| Ok(r))?;
    let beta_var = cs.new_input_variable(|| Ok(beta))?;

    cs.enforce_constraint(lc!() + beta_var, lc!() + t1_x_var, lc!() + r_var - t2_x_var)?;

    Ok((r_var, r))
}

/// R_T-eVRF-DH with the paper's §5.3 Weierstrass gadget (3 constraints per bit).
///
/// All points are passed as Edwards points and converted to short Weierstrass
/// coordinates; S.X, T1.X, T2.X are Weierstrass x-coordinates, as in the paper.
/// The public bases (g_in, PK2, H1, H2) and the correction generator h_c enter
/// only through precomputed tables, i.e. as matrix coefficients.
/// Public inputs, in allocation order: [1, PK1.X, PK1.Y, beta].
///
/// Written with AI assistance (Claude), following Golden §5.3.
pub fn evrf_dh_sw(
    cs: ConstraintSystemRef<ark_ed_on_bls12_381_bandersnatch::Fq>,
    g_in: &ark_ed_on_bls12_381_bandersnatch::EdwardsAffine,
    pk1: &ark_ed_on_bls12_381_bandersnatch::EdwardsAffine,
    pk2: &ark_ed_on_bls12_381_bandersnatch::EdwardsAffine,
    h1: &ark_ed_on_bls12_381_bandersnatch::EdwardsAffine,
    h2: &ark_ed_on_bls12_381_bandersnatch::EdwardsAffine,
    h_c: &ark_ed_on_bls12_381_bandersnatch::EdwardsAffine,
    beta: ark_ed_on_bls12_381_bandersnatch::Fq,
    sk: ark_ed_on_bls12_381_bandersnatch::Fq,
) -> Result<(Variable, ark_ed_on_bls12_381_bandersnatch::Fq), SynthesisError> {
    use crate::sw::{point_to_sw, scalar_mul_sw_fixed, sw_fixed_table};

    // (0) PK1 = g_in^sk
    let sk_var = cs.new_witness_variable(|| Ok(sk))?;
    let (sk_bits_var, sk_bits) = bit_decompose(cs.clone(), sk, sk_var)?;
    let nbits = sk_bits_var.len();

    let t_g = sw_fixed_table(g_in, h_c, nbits);
    let (px_var, py_var, _, _) = scalar_mul_sw_fixed(cs.clone(), &sk_bits_var, &sk_bits, &t_g)?;

    let (pk1_x, pk1_y) = point_to_sw(pk1);
    let pk1_x_var = cs.new_input_variable(|| Ok(pk1_x))?;
    let pk1_y_var = cs.new_input_variable(|| Ok(pk1_y))?;
    cs.enforce_constraint(lc!() + pk1_x_var, lc!() + Variable::One, lc!() + px_var)?;
    cs.enforce_constraint(lc!() + pk1_y_var, lc!() + Variable::One, lc!() + py_var)?;

    // (1) S = PK2^sk ; (2)-(3) k = S.X
    let t_pk2 = sw_fixed_table(pk2, h_c, nbits);
    let (sx_var, _, sx, _) = scalar_mul_sw_fixed(cs.clone(), &sk_bits_var, &sk_bits, &t_pk2)?;
    let (k_bits_var, k_bits) = bit_decompose(cs.clone(), sx, sx_var)?;

    // (4)-(5) T_j = H_j^k ; (6)-(7) r_j = T_j.X
    let t_h1 = sw_fixed_table(h1, h_c, nbits);
    let t_h2 = sw_fixed_table(h2, h_c, nbits);
    let (t1x_var, _, t1x, _) = scalar_mul_sw_fixed(cs.clone(), &k_bits_var, &k_bits, &t_h1)?;
    let (t2x_var, _, t2x, _) = scalar_mul_sw_fixed(cs.clone(), &k_bits_var, &k_bits, &t_h2)?;

    // (8) r = beta * r1 + r2
    let r = beta * t1x + t2x;
    let r_var = cs.new_witness_variable(|| Ok(r))?;
    let beta_var = cs.new_input_variable(|| Ok(beta))?;
    cs.enforce_constraint(lc!() + beta_var, lc!() + t1x_var, lc!() + r_var - t2x_var)?;

    Ok((r_var, r))
}
