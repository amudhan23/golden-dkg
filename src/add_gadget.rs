use ark_ff::{BigInteger, PrimeField};
use ark_relations::{
    lc,
    r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError, Variable},
};

pub fn edwards_add<F: PrimeField>(
    cs: ConstraintSystemRef<F>,
    x1: Variable,
    y1: Variable,
    x2: Variable,
    y2: Variable,
    x1_val: F,
    y1_val: F,
    x2_val: F,
    y2_val: F,
    a: F, // curve parameter
    d: F, // curve parameter
) -> Result<(Variable, Variable, F, F), SynthesisError> {
    let x1_x2_val = x1_val * x2_val;
    let x1_x2 = cs.new_witness_variable(|| Ok(x1_x2_val))?;

    cs.enforce_constraint(lc!() + x1, lc!() + x2, lc!() + x1_x2)?;

    let y1_y2_val = y1_val * y2_val;
    let y1_y2 = cs.new_witness_variable(|| Ok(y1_y2_val))?;

    cs.enforce_constraint(lc!() + y1, lc!() + y2, lc!() + y1_y2)?;

    let x1_x2_y1_y2_val = x1_x2_val * y1_y2_val;

    let x1_x2_y1_y2 = cs.new_witness_variable(|| Ok(x1_x2_y1_y2_val))?;

    cs.enforce_constraint(lc!() + x1_x2, lc!() + y1_y2, lc!() + x1_x2_y1_y2)?;

    let x1_y2_val = x1_val * y2_val;
    let x1_y2 = cs.new_witness_variable(|| Ok(x1_y2_val))?;
    cs.enforce_constraint(lc!() + x1, lc!() + y2, lc!() + x1_y2)?;

    let y1_x2_val = y1_val * x2_val;

    let y1_x2 = cs.new_witness_variable(|| Ok(y1_x2_val))?;
    cs.enforce_constraint(lc!() + y1, lc!() + x2, lc!() + y1_x2)?;

    let denom_x = F::one() + d * x1_x2_y1_y2_val;
    let x3_val = (x1_y2_val + y1_x2_val) * denom_x.inverse().unwrap();
    let x3 = cs.new_witness_variable(|| Ok(x3_val))?;

    let denom_y = F::one() - d * x1_x2_y1_y2_val;
    let y3_val = (y1_y2_val - a * x1_x2_val) * denom_y.inverse().unwrap();
    let y3 = cs.new_witness_variable(|| Ok(y3_val))?;

    cs.enforce_constraint(
        lc!() + x3,
        lc!() + Variable::One + (d, x1_x2_y1_y2),
        lc!() + x1_y2 + y1_x2,
    )?;
    cs.enforce_constraint(
        lc!() + y3,
        lc!() + Variable::One - (d, x1_x2_y1_y2),
        lc!() + y1_y2 - (a, x1_x2),
    )?;

    Ok((x3, y3, x3_val, y3_val))
}

pub fn edwards_add_newConstraint<F: PrimeField>(
    cs: ConstraintSystemRef<F>,
    x1: Variable,
    y1: Variable,
    x1_val: F,
    y1_val: F,
    x2_val: F,
    y2_val: F,
    bit: F,
    bit_var: Variable,
    a: F, // curve parameter
    d: F, // curve parameter
) -> Result<(Variable, Variable, F, F), SynthesisError> {
    let m1 = bit * x1_val;
    let m1_var = cs.new_witness_variable(|| Ok(m1))?;

    cs.enforce_constraint(lc!() + bit_var, lc!() + x1, lc!() + m1_var)?;

    let m2 = bit * y1_val;
    let m2_var = cs.new_witness_variable(|| Ok(m2))?;

    cs.enforce_constraint(lc!() + bit_var, lc!() + y1, lc!() + m2_var)?;

    let p = (x2_val * m1) * (y1_val + (y2_val - F::one()) * m2);
    let p_var = cs.new_witness_variable(|| Ok(p))?;

    cs.enforce_constraint(
        lc!() + (x2_val, m1_var),
        lc!() + y1 + (y2_val - F::one(), m2_var),
        lc!() + p_var,
    )?;

    let denom_x = F::one() + d * p;
    let x3_val = (x1_val + (y2_val - F::one()) * m1 + x2_val * m2) * denom_x.inverse().unwrap();
    let x3 = cs.new_witness_variable(|| Ok(x3_val))?;

    let denom_y = F::one() - d * p;
    let y3_val = (y1_val + (y2_val - F::one()) * m2 - a * x2_val * m1) * denom_y.inverse().unwrap();
    let y3 = cs.new_witness_variable(|| Ok(y3_val))?;

    cs.enforce_constraint(
        lc!() + x3,
        lc!() + Variable::One + (d, p_var),
        lc!() + x1 + (y2_val - F::one(), m1_var) + (x2_val, m2_var),
    )?;
    cs.enforce_constraint(
        lc!() + y3,
        lc!() + Variable::One - (d, p_var),
        lc!() + y1 + (y2_val - F::one(), m2_var) - (a * x2_val, m1_var),
    )?;

    Ok((x3, y3, x3_val, y3_val))
}

pub fn edwards_add_native<F: PrimeField>(x1: F, y1: F, x2: F, y2: F, a: F, d: F) -> (F, F) {
    let t = d * x1 * x2 * y1 * y2;
    let x3 = (x1 * y2 + y1 * x2) * (F::one() + t).inverse().unwrap();
    let y3 = (y1 * y2 - a * x1 * x2) * (F::one() - t).inverse().unwrap();
    (x3, y3)
}
