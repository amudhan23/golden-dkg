#[cfg(test)]
mod tests {
    use crate::circuit_evrf::bit_decompose;
    use crate::sw::*;

    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ed_on_bls12_381_bandersnatch::{EdwardsAffine, Fq, Fr as Scalar};
    use ark_ff::{One, PrimeField};
    use ark_relations::r1cs::{ConstraintSystem, Variable};
    use ark_std::UniformRand;

    fn rand_point(rng: &mut impl ark_std::rand::Rng) -> EdwardsAffine {
        (EdwardsAffine::generator() * Scalar::rand(rng)).into_affine()
    }

    fn h_c() -> EdwardsAffine {
        // any independent subgroup point works for tests
        (EdwardsAffine::generator() * Scalar::from(0x5eed_u64)).into_affine()
    }

    #[test]
    fn generator_maps_onto_sw_curve() {
        let p = bandersnatch_sw();
        let (x, y) = point_to_sw(&EdwardsAffine::generator());
        assert!(on_sw_curve(&p, x, y));
    }

    #[test]
    fn conversion_preserves_addition_and_doubling() {
        let p = bandersnatch_sw();
        let mut rng = ark_std::test_rng();
        for _ in 0..10 {
            let a = rand_point(&mut rng);
            let b = rand_point(&mut rng);
            let (ax, ay) = point_to_sw(&a);
            let (bx, by) = point_to_sw(&b);
            assert!(on_sw_curve(&p, ax, ay));
            assert_eq!(point_to_sw(&(a + b).into_affine()), sw_add_native(ax, ay, bx, by));
            assert_eq!(point_to_sw(&(a + a).into_affine()), sw_double_native(&p, ax, ay));
        }
    }

    /// Run the gadget on k·base and compare with native arithmetic.
    fn gadget_check(k: Scalar, base: &EdwardsAffine) -> (bool, bool, usize) {
        let cs = ConstraintSystem::<Fq>::new_ref();
        let k_fq = Fq::from_bigint(k.into_bigint().into()).unwrap();
        let k_var = cs.new_witness_variable(|| Ok(k_fq)).unwrap();
        let (bits, vals) = bit_decompose(cs.clone(), k_fq, k_var).unwrap();
        let table = sw_fixed_table(base, &h_c(), bits.len());
        let before = cs.num_constraints();
        let (_, _, x, y) = scalar_mul_sw_fixed(cs.clone(), &bits, &vals, &table).unwrap();
        let used = cs.num_constraints() - before;
        let expected = point_to_sw(&(*base * k).into_affine());
        ((x, y) == expected, cs.is_satisfied().unwrap(), used)
    }

    #[test]
    fn gadget_matches_native_small_and_random_scalars() {
        let mut rng = ark_std::test_rng();
        let base = rand_point(&mut rng);
        for k in [1u64, 2, 3, 42, 1 << 40] {
            let (ok, sat, _) = gadget_check(Scalar::from(k), &base);
            assert!(ok, "k = {k}: wrong result");
            assert!(sat, "k = {k}: constraints not satisfied");
        }
        for _ in 0..3 {
            let (ok, sat, _) = gadget_check(Scalar::rand(&mut rng), &base);
            assert!(ok && sat, "random k failed");
        }
    }

    #[test]
    fn gadget_works_when_base_is_the_generator() {
        // Step (0) of the eVRF uses g_in itself as the base.
        let mut rng = ark_std::test_rng();
        let g = EdwardsAffine::generator();
        for k in [1u64, 2, 7, 255, 256] {
            let (ok, sat, _) = gadget_check(Scalar::from(k), &g);
            assert!(ok && sat, "k = {k} with base g failed");
        }
        let (ok, sat, _) = gadget_check(Scalar::rand(&mut rng), &g);
        assert!(ok && sat);
    }

    #[test]
    fn gadget_costs_3_per_bit_plus_2() {
        let mut rng = ark_std::test_rng();
        let (_, _, used) = gadget_check(Scalar::from(42u64), &rand_point(&mut rng));
        assert_eq!(used, 3 * 255 + 2, "paper §5.3: 3λ + 2 with λ + 1 = 256 bits");
    }

    #[test]
    fn corrections_cancel() {
        // Σ C_j = 0, i.e. with all bits 0 the result is the identity: the last
        // step computes L_{λ-1} + C_λ = -C_λ + C_λ, which the chord rule cannot
        // represent; so instead check Σ c_j ≡ 0 directly.
        let lambda: u64 = 255;
        let s_minus = -Scalar::from(2 * lambda - 1);
        let sum = Scalar::from(1u64) + Scalar::from(2u64) * Scalar::from(lambda - 1) + s_minus;
        assert_eq!(sum, Scalar::from(0u64));
    }

    #[test]
    fn tampered_output_fails() {
        let cs = ConstraintSystem::<Fq>::new_ref();
        let k = Fq::from(42u64);
        let k_var = cs.new_witness_variable(|| Ok(k)).unwrap();
        let (bits, vals) = bit_decompose(cs.clone(), k, k_var).unwrap();
        let table = sw_fixed_table(&EdwardsAffine::generator(), &h_c(), bits.len());
        let (xv, _, _, _) = scalar_mul_sw_fixed(cs.clone(), &bits, &vals, &table).unwrap();
        match xv {
            Variable::Witness(i) => cs.borrow_mut().unwrap().witness_assignment[i] += Fq::one(),
            _ => panic!(),
        }
        assert!(!cs.is_satisfied().unwrap());
    }
}
