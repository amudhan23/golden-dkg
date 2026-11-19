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
