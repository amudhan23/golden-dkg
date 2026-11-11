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

        let bits = bit_decompose(cs.clone(), scalar).unwrap();

        assert_eq!(bits.len(), 256);
        assert!(
            cs.is_satisfied().unwrap(),
            "bit decompostion constraints not satisfied"
        );
        println!("constraints: {}", cs.num_constraints());
    }
}
