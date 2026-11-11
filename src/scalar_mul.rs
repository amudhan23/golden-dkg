use ark_ff::{BigInteger, PrimeField};
use ark_relations::{
    lc,
    r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError, Variable},
};

use crate::add_gadget::edwards_add;

pub fn scalar_mul<F: PrimeField>(
    cs: ConstraintSystemRef<F>,
    bits: &[Variable], // from bit_decompose, LSB first
    bit_vals: &[F],    // the actual bit values
    base_x: F,
    base_y: F, // the base point P (public, known)
    a: F,
    d: F, // curve params
) -> Result<(Variable, Variable, F, F), SynthesisError> {
    let mut acc_x_val = F::zero();
    let mut acc_y_val = F::one();
    let mut acc_x_var = cs.new_witness_variable(|| Ok(acc_x_val))?;
    let mut acc_y_var = cs.new_witness_variable(|| Ok(acc_y_val))?;

    let p_var_x = cs.new_witness_variable(|| Ok(base_x))?;
    let p_var_y = cs.new_witness_variable(|| Ok(base_y))?;

    for i in (0..bits.len()).rev() {
        let (dbl_x_var, dbl_y_var, dbl_x_val, dbl_y_val) = edwards_add(
            cs.clone(),
            acc_x_var,
            acc_y_var,
            acc_x_var,
            acc_y_var,
            acc_x_val,
            acc_y_val,
            acc_x_val,
            acc_y_val,
            a,
            d,
        )?;

        let (add_x_var, add_y_var, add_x_val, add_y_val) = edwards_add(
            cs.clone(),
            dbl_x_var,
            dbl_y_var,
            p_var_x,
            p_var_y,
            dbl_x_val,
            dbl_y_val,
            base_x,
            base_y,
            a,
            d,
        )?;

        (acc_x_var, acc_y_var, acc_x_val, acc_y_val) = conditional_select(
            cs.clone(),
            bits[i],
            bit_vals[i],
            add_x_var,
            add_y_var,
            add_x_val,
            add_y_val,
            dbl_x_var,
            dbl_y_var,
            dbl_x_val,
            dbl_y_val,
        )?;
    }

    Ok((acc_x_var, acc_y_var, acc_x_val, acc_y_val))
}

pub fn conditional_select<F: PrimeField>(
    cs: ConstraintSystemRef<F>,
    bit: Variable, // the selector bit (0 or 1)
    bit_val: F,
    // "if 1" point
    x_t: Variable,
    y_t: Variable,
    x_t_val: F,
    y_t_val: F,
    // "if 0" point
    x_a: Variable,
    y_a: Variable,
    x_a_val: F,
    y_a_val: F,
) -> Result<(Variable, Variable, F, F), SynthesisError> {
    let out_x_val = bit_val * x_t_val + (F::one() - bit_val) * x_a_val;
    let out_y_val = bit_val * y_t_val + (F::one() - bit_val) * y_a_val;

    let out_x_var = cs.new_witness_variable(|| Ok(out_x_val))?;
    let out_y_var = cs.new_witness_variable(|| Ok(out_y_val))?;

    cs.enforce_constraint(lc!() + bit, lc!() + x_t - x_a, lc!() + out_x_var - x_a)?;

    cs.enforce_constraint(lc!() + bit, lc!() + y_t - y_a, lc!() + out_y_var - y_a)?;

    Ok((out_x_var, out_y_var, out_x_val, out_y_val))
}
