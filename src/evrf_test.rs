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
        // fixed-base circuit: ~5.6k constraints + ~5.6k witness variables -> 16,384 gates
        assert_eq!(p1.bp.n, 16_384);
        assert_eq!(p1.bp.g.len(), p1.bp.n);
    }

    /// Alice -> Bob and Bob -> Alice must produce the SAME pad r (DH symmetry +
    /// sorted hash inputs), R must equal r*U, and Alice's proofs must verify
    /// against her R via `verify`, and wrong R / swapped keys / wrong msg /
    /// non-subgroup key are rejected. Slow (~9 min):
    ///     cargo test --release evaluate_both_sides -- --ignored --nocapture
    #[test]
    #[ignore]
    fn evaluate_both_sides_agree_and_proofs_verify() {
        use crate::evrf::evaluate;
        use ark_ed_on_bls12_381_bandersnatch::Fq;
        use std::time::Instant;

        let params = setup_evrf();
        let mut rng = ark_std::test_rng();
        let g = EdwardsAffine::generator();
        let sk_a = Scalar::rand(&mut rng);
        let sk_b = Scalar::rand(&mut rng);
        let pk_a = (g * sk_a).into_affine();
        let pk_b = (g * sk_b).into_affine();
        let msg = b"session-1";

        let t = Instant::now();
        let (r_ab, big_r_ab, proof_ab) = evaluate(&params, sk_a, &pk_a, &pk_b, msg, &mut rng);
        println!("evaluate (Alice -> Bob): {:.1} s", t.elapsed().as_secs_f64());
        let (r_ba, big_r_ba, _) = evaluate(&params, sk_b, &pk_b, &pk_a, msg, &mut rng);

        assert_eq!(r_ab, r_ba, "both parties must derive the same pad r");
        assert_eq!(big_r_ab, big_r_ba);
        assert_eq!(big_r_ab, params.bp.u * r_ab, "R must be r * U");

        use crate::evrf::verify;

        let t = Instant::now();
        let ok = verify(&params, &pk_a, &pk_b, msg, big_r_ab, &proof_ab);
        println!("verify: {:.1} s", t.elapsed().as_secs_f64());
        assert!(ok, "Alice's eVRF proof must verify against her R");

        // ── Rejections ──
        // The fake-R attack: an R for a different pad.
        let fake_r = params.bp.u * (r_ab + Fq::from(1u64));
        assert!(!verify(&params, &pk_a, &pk_b, msg, fake_r, &proof_ab), "wrong R must fail");

        // Roles swapped: the proof says Alice is the prover (PK1 = pk_a).
        assert!(!verify(&params, &pk_b, &pk_a, msg, big_r_ab, &proof_ab), "swapped keys must fail");

        // Different session message -> different H1, H2.
        assert!(!verify(&params, &pk_a, &pk_b, b"session-2", big_r_ab, &proof_ab), "wrong msg must fail");

        // A key outside the prime-order subgroup: (0, -1) is on the curve with order 2.
        let bad = EdwardsAffine::new_unchecked(Fq::from(0u64), -Fq::from(1u64));
        assert!(bad.is_on_curve() && !bad.is_in_correct_subgroup_assuming_on_curve());
        assert!(!verify(&params, &bad, &pk_b, msg, big_r_ab, &proof_ab), "non-subgroup key must fail");
    }
}
