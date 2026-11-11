use ark_ff::{BigInteger, PrimeField};
use ark_relations::{
    lc,
    r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError, Variable},
};

pub fn bit_decompose<F: PrimeField>(
    cs: ConstraintSystemRef<F>,
    scalar: F,
) -> Result<Vec<Variable>, SynthesisError> {
    let bigint = scalar.into_bigint();
    let mut big_vec = Vec::with_capacity(256);

    for i in 0..256 {
        let bit_val = if bigint.get_bit(i) {
            F::one()
        } else {
            F::zero()
        };

        let bit_var = cs.new_witness_variable(|| Ok(bit_val))?;

        big_vec.push(bit_var);
    }

    for i in 0..256 {
        cs.enforce_constraint(
            lc!() + big_vec[i],
            lc!() + big_vec[i] - Variable::One,
            lc!(),
        )?;
    }

    let mut reconstruction = lc!();
    let mut coeff = F::one();
    let two = F::from(2u64);

    for i in 0..256 {
        reconstruction = reconstruction + (coeff, big_vec[i]);
        coeff *= two;
    }

    let k_var = cs.new_witness_variable(|| Ok(scalar))?;

    cs.enforce_constraint(reconstruction, lc!() + Variable::One, lc!() + k_var)?;

    Ok(big_vec)
}
