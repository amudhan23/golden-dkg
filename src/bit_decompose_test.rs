#[cfg(test)]
mod tests {
    use crate::circuit_evrf::bit_decompose;
    use ark_bls12_381::Fr;
    use ark_ff::PrimeField;
    use ark_relations::r1cs::ConstraintSystem;

    #[test]
    fn bit_decompose_satisfies() {
        let cs = ConstraintSystem::<Fr>::new_ref();
        let scalar = Fr::from(7u64);

        let scalar_var = cs.new_witness_variable(|| Ok(scalar)).unwrap();
        let (bits, bit_vals) = bit_decompose(cs.clone(), scalar, scalar_var).unwrap();
        assert_eq!(bit_vals.len(), 256);

        assert_eq!(bits.len(), 256);
        assert!(
            cs.is_satisfied().unwrap(),
            "bit decompostion constraints not satisfied"
        );
        println!("constraints: {}", cs.num_constraints());
    }

    #[test]
    fn bit_decompose_rejects_mismatched_scalar_var() {
        // Bits are computed from 7, but the linked variable holds 8.
        // The reconstruction constraint must catch this.
        let cs = ConstraintSystem::<Fr>::new_ref();
        let scalar = Fr::from(7u64);
        let wrong_var = cs.new_witness_variable(|| Ok(Fr::from(8u64))).unwrap();

        bit_decompose(cs.clone(), scalar, wrong_var).unwrap();

        assert!(
            !cs.is_satisfied().unwrap(),
            "mismatched scalar_var should not satisfy the constraints"
        );
    }
}
