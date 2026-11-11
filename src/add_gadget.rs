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
