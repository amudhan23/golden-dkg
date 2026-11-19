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

    /// Alice -> Bob and Bob -> Alice must produce the SAME pad r (DH symmetry +
    /// sorted hash inputs), R must equal r*U, and Alice's proofs must verify
    /// against her R. Slow (two proofs, ~4 min):
    ///     cargo test --release evaluate_both_sides -- --ignored --nocapture
    #[test]
    #[ignore]
    fn evaluate_both_sides_agree_and_proofs_verify() {
        use crate::evrf::evaluate;
        use crate::schnorr;
        use ark_ed_on_bls12_381_bandersnatch::Fq;
        use ark_relations::r1cs::{ConstraintSystem, SynthesisMode, Variable};
        use merlin::Transcript;
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

        // Verifier side (a preview of `verify`): own matrices in Setup mode.
        let (h1, h2) = hash_points(msg, &pk_a, &pk_b);
        let cs = ConstraintSystem::<Fq>::new_ref();
        cs.set_mode(SynthesisMode::Setup);
        let (r_var, _) = crate::evrf_dh_circuit::evrf_dh(
            cs.clone(),
            params.g_in.x, params.g_in.y,
            pk_a.x, pk_a.y,
            pk_b.x, pk_b.y,
            h1.x, h1.y,
            h2.x, h2.y,
            params.beta,
            Fq::from(0u64),
            params.a, params.d,
        )
        .unwrap();
        cs.finalize();
        let mats = cs.to_matrices().unwrap();
        let r_index = match r_var { Variable::Witness(i) => i, _ => panic!() };
        let instance = vec![
            Fq::from(1u64),
            params.g_in.x, params.g_in.y,
            pk_a.x, pk_a.y,
            pk_b.x, pk_b.y,
            h1.x, h1.y,
            h2.x, h2.y,
            params.beta,
        ];

        let t = Instant::now();
        let mut tr = Transcript::new(b"unit-golden");
        let bp_ok = bulletproofs_ipa::protocol::verify(
            &mats, &instance, &params.bp, &[r_index], &[big_r_ab], &proof_ab.bp, &mut tr,
        );
        let pok_ok = schnorr::verify(&proof_ab.pok, params.bp.u, big_r_ab, &mut tr);
        println!("verify (Bulletproof + Schnorr): {:.1} s", t.elapsed().as_secs_f64());
        assert!(bp_ok, "Alice's Bulletproof must verify against her R");
        assert!(pok_ok, "Alice's Schnorr proof for R must verify");
    }
}
