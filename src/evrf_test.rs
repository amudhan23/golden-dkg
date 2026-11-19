#[cfg(test)]
mod tests {
    use crate::evrf::{hash_points, setup_evrf};

    use ark_ec::twisted_edwards::TECurveConfig;
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ed_on_bls12_381_bandersnatch::{BandersnatchConfig, EdwardsAffine, Fr as Scalar};
    use ark_std::UniformRand;

    fn two_keys() -> (EdwardsAffine, EdwardsAffine) {
        let mut rng = ark_std::test_rng();
        let g = EdwardsAffine::generator();
        let a = (g * Scalar::rand(&mut rng)).into_affine();
        let b = (g * Scalar::rand(&mut rng)).into_affine();
        assert_ne!(a, b);
        (a, b)
    }

    #[test]
    fn hash_points_is_symmetric() {
        // Alice and Bob must derive the same H1, H2 whoever calls.
        let (pk_a, pk_b) = two_keys();
        assert_eq!(hash_points(b"msg", &pk_a, &pk_b), hash_points(b"msg", &pk_b, &pk_a));
    }

    #[test]
    fn hash_points_h1_and_h2_differ() {
        let (pk_a, pk_b) = two_keys();
        let (h1, h2) = hash_points(b"msg", &pk_a, &pk_b);
        assert_ne!(h1, h2, "different domain tags must give independent points");
    }

    #[test]
    fn hash_points_depend_on_inputs() {
        let (pk_a, pk_b) = two_keys();
        let base = hash_points(b"msg", &pk_a, &pk_b);
        assert_ne!(base, hash_points(b"other msg", &pk_a, &pk_b), "msg must matter");
        let pk_c = (pk_b + EdwardsAffine::generator()).into_affine();
        assert_ne!(base, hash_points(b"msg", &pk_a, &pk_c), "keys must matter");
    }

    #[test]
    fn hash_points_are_in_prime_order_subgroup() {
        // Required for the Edwards addition formula to be exception-free.
        let (pk_a, pk_b) = two_keys();
        let (h1, h2) = hash_points(b"msg", &pk_a, &pk_b);
        for h in [h1, h2] {
            assert!(h.is_on_curve());
            assert!(h.is_in_correct_subgroup_assuming_on_curve());
            assert!(!h.is_zero());
        }
    }

    #[test]
    fn setup_is_deterministic_and_correct() {
        let p1 = setup_evrf();
        let p2 = setup_evrf();
        assert_eq!(p1.beta, p2.beta, "beta must be the same for everyone");
        assert_eq!(p1.bp.g, p2.bp.g, "generators must be the same for everyone");
        assert_eq!(p1.bp.u, p2.bp.u);
        assert_eq!(p1.g_in, EdwardsAffine::generator());
        assert_eq!(p1.a, BandersnatchConfig::COEFF_A);
        assert_eq!(p1.d, BandersnatchConfig::COEFF_D);
        // 16,909 constraints + 16,906 witness variables -> 65,536 gates
        assert_eq!(p1.bp.n, 65_536);
        assert_eq!(p1.bp.g.len(), p1.bp.n);
    }
}
