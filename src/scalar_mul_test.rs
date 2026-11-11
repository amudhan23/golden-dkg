#[cfg(test)]
mod tests {
    use crate::circuit_evrf::bit_decompose;
    use crate::scalar_mul::{conditional_select, scalar_mul};

    use ark_ec::twisted_edwards::TECurveConfig;
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ed_on_bls12_381_bandersnatch::{BandersnatchConfig, EdwardsAffine, Fq};
    use ark_ff::{BigInteger, PrimeField, One, Zero};
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::UniformRand;

    /// Helper: extract bit values (as field elements) from a scalar, LSB first
    fn bit_vals_from_scalar(scalar: Fq) -> Vec<Fq> {
        let bigint = scalar.into_bigint();
        (0..256)
            .map(|i| {
                if bigint.get_bit(i) {
                    Fq::one()
                } else {
                    Fq::zero()
                }
            })
            .collect()
    }

    #[test]
    fn conditional_select_picks_first_when_bit_is_one() {
        let cs = ConstraintSystem::<Fq>::new_ref();

        let x_t_val = Fq::from(10u64);
        let y_t_val = Fq::from(20u64);
        let x_a_val = Fq::from(30u64);
        let y_a_val = Fq::from(40u64);

        let bit = cs.new_witness_variable(|| Ok(Fq::one())).unwrap();
        let x_t = cs.new_witness_variable(|| Ok(x_t_val)).unwrap();
        let y_t = cs.new_witness_variable(|| Ok(y_t_val)).unwrap();
        let x_a = cs.new_witness_variable(|| Ok(x_a_val)).unwrap();
        let y_a = cs.new_witness_variable(|| Ok(y_a_val)).unwrap();

        let (_, _, out_x_val, out_y_val) = conditional_select(
            cs.clone(),
            bit,
            Fq::one(),
            x_t, y_t, x_t_val, y_t_val,
            x_a, y_a, x_a_val, y_a_val,
        )
        .unwrap();

        assert_eq!(out_x_val, x_t_val, "should select x_t when bit=1");
        assert_eq!(out_y_val, y_t_val, "should select y_t when bit=1");
        assert!(cs.is_satisfied().unwrap(), "constraints not satisfied for bit=1");
        assert_eq!(cs.num_constraints(), 2, "conditional_select should have 2 constraints");
    }

    #[test]
    fn conditional_select_picks_second_when_bit_is_zero() {
        let cs = ConstraintSystem::<Fq>::new_ref();

        let x_t_val = Fq::from(10u64);
        let y_t_val = Fq::from(20u64);
        let x_a_val = Fq::from(30u64);
        let y_a_val = Fq::from(40u64);

        let bit = cs.new_witness_variable(|| Ok(Fq::zero())).unwrap();
        let x_t = cs.new_witness_variable(|| Ok(x_t_val)).unwrap();
        let y_t = cs.new_witness_variable(|| Ok(y_t_val)).unwrap();
        let x_a = cs.new_witness_variable(|| Ok(x_a_val)).unwrap();
        let y_a = cs.new_witness_variable(|| Ok(y_a_val)).unwrap();

        let (_, _, out_x_val, out_y_val) = conditional_select(
            cs.clone(),
            bit,
            Fq::zero(),
            x_t, y_t, x_t_val, y_t_val,
            x_a, y_a, x_a_val, y_a_val,
        )
        .unwrap();

        assert_eq!(out_x_val, x_a_val, "should select x_a when bit=0");
        assert_eq!(out_y_val, y_a_val, "should select y_a when bit=0");
        assert!(cs.is_satisfied().unwrap(), "constraints not satisfied for bit=0");
    }

    #[test]
    fn scalar_mul_known_scalar() {
        // Test: compute k * G using scalar_mul gadget and compare with arkworks native
        let cs = ConstraintSystem::<Fq>::new_ref();

        // Bandersnatch curve params
        let a: Fq = BandersnatchConfig::COEFF_A;
        let d: Fq = BandersnatchConfig::COEFF_D;

        // Use a small known scalar
        let k = Fq::from(42u64);

        // Get the generator point in affine twisted Edwards form
        let generator: EdwardsAffine = BandersnatchConfig::GENERATOR;
        let base_x: Fq = generator.x;
        let base_y: Fq = generator.y;

        // Compute expected result natively: k * G
        // We need to convert k (Fq) to the scalar field (Fr) for native mul
        // But Bandersnatch's base field IS BLS12-381's Fr, and its scalar field
        // is the Bandersnatch subgroup order. For native EC mul, arkworks uses
        // the scalar field. We need to be careful here.
        //
        // Actually: bit_decompose decomposes k as an Fq element (256 bits).
        // scalar_mul treats those bits as the scalar for double-and-add.
        // The native equivalent is interpreting k's bit pattern as an integer
        // and doing that many additions of G.
        // For small k (like 42), this is just 42*G in the group.
        let expected = EdwardsAffine::generator() * ark_ed_on_bls12_381_bandersnatch::Fr::from(42u64);
        let expected_affine: EdwardsAffine = expected.into_affine();

        // Bit decompose k
        let bit_values = bit_vals_from_scalar(k);
        let bits = bit_decompose(cs.clone(), k).unwrap();

        // Run scalar_mul
        let (_, _, result_x, result_y) =
            scalar_mul(cs.clone(), &bits, &bit_values, base_x, base_y, a, d).unwrap();

        // Check the result matches native computation
        assert_eq!(
            result_x, expected_affine.x,
            "scalar_mul x-coordinate doesn't match native"
        );
        assert_eq!(
            result_y, expected_affine.y,
            "scalar_mul y-coordinate doesn't match native"
        );

        // Check all constraints are satisfied
        assert!(
            cs.is_satisfied().unwrap(),
            "scalar_mul constraints not satisfied"
        );

        let num_constraints = cs.num_constraints();
        println!("scalar_mul constraints: {}", num_constraints);
        // Expected: 257 (bit_decompose) + 256 * 16 (scalar_mul) = 4353
        println!(
            "  bit_decompose: 257, scalar_mul loop: {} * 16 = {}",
            bits.len(),
            bits.len() * 16
        );
    }

    #[test]
    fn scalar_mul_scalar_one() {
        // k=1 should give G itself
        let cs = ConstraintSystem::<Fq>::new_ref();

        let a: Fq = BandersnatchConfig::COEFF_A;
        let d: Fq = BandersnatchConfig::COEFF_D;
        let k = Fq::from(1u64);

        let generator: EdwardsAffine = BandersnatchConfig::GENERATOR;
        let base_x: Fq = generator.x;
        let base_y: Fq = generator.y;

        let bit_values = bit_vals_from_scalar(k);
        let bits = bit_decompose(cs.clone(), k).unwrap();

        let (_, _, result_x, result_y) =
            scalar_mul(cs.clone(), &bits, &bit_values, base_x, base_y, a, d).unwrap();

        // 1 * G = G
        assert_eq!(result_x, base_x, "1*G should give G.x");
        assert_eq!(result_y, base_y, "1*G should give G.y");
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn scalar_mul_random_scalar() {
        // Test with a random scalar
        let cs = ConstraintSystem::<Fq>::new_ref();

        let a: Fq = BandersnatchConfig::COEFF_A;
        let d: Fq = BandersnatchConfig::COEFF_D;

        let mut rng = ark_std::test_rng();
        let k_scalar = ark_ed_on_bls12_381_bandersnatch::Fr::rand(&mut rng);

        // Convert scalar field element to base field for bit decomposition
        // We use the scalar's BigInt representation directly
        let k_bigint = k_scalar.into_bigint();
        let k_fq = Fq::from_bigint(k_bigint.into()).unwrap();

        let generator: EdwardsAffine = BandersnatchConfig::GENERATOR;
        let base_x: Fq = generator.x;
        let base_y: Fq = generator.y;

        // Native computation
        let expected = EdwardsAffine::generator() * k_scalar;
        let expected_affine: EdwardsAffine = expected.into_affine();

        // Circuit computation
        let bit_values = bit_vals_from_scalar(k_fq);
        let bits = bit_decompose(cs.clone(), k_fq).unwrap();

        let (_, _, result_x, result_y) =
            scalar_mul(cs.clone(), &bits, &bit_values, base_x, base_y, a, d).unwrap();

        assert_eq!(result_x, expected_affine.x, "random scalar_mul x mismatch");
        assert_eq!(result_y, expected_affine.y, "random scalar_mul y mismatch");
        assert!(cs.is_satisfied().unwrap());

        println!(
            "random scalar_mul: {} constraints",
            cs.num_constraints()
        );
    }
}
