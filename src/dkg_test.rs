#[cfg(test)]
mod tests {
    use crate::dkg::round0;
    use crate::evrf::{setup_evrf, verify, EvrfParams};
    use crate::feldman::{compute_r, share_commitment, verify_share};

    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ed_on_bls12_381_bandersnatch::{EdwardsAffine, Fr as Scalar};
    use ark_std::UniformRand;

    fn parties(n: usize) -> (Vec<Scalar>, Vec<EdwardsAffine>) {
        let mut rng = ark_std::test_rng();
        let g = EdwardsAffine::generator();
        let sks: Vec<Scalar> = (0..n).map(|_| Scalar::rand(&mut rng)).collect();
        let pks = sks.iter().map(|sk| (g * sk).into_affine()).collect();
        (sks, pks)
    }

    /// One dealing from party 1 in a 3-party, threshold-2 setup:
    /// structure, public ciphertext check, eVRF proof, and decryption.
    /// Slow (2 proofs + 1 verify, ~30 s):
    ///     cargo test --release round0_dealing -- --ignored --nocapture
    #[test]
    #[ignore]
    fn round0_dealing_is_well_formed_and_decryptable() {
        let params: EvrfParams = setup_evrf();
        let (n, t) = (3usize, 2usize);
        let (sks, pks) = parties(n);
        let mut rng = ark_std::test_rng();
        let g = params.bp.u;

        let (dealing, own_share) = round0(&params, n, t, 1, sks[0], &pks, &mut rng);

        assert_eq!(dealing.dealer, 1);
        assert_eq!(dealing.commitments.len(), t, "t Feldman commitments");
        assert_eq!(dealing.entries.len(), n - 1, "one entry per other party");
        assert!(dealing.entries.iter().all(|e| e.recipient != 1), "no entry for the dealer");
        assert!(verify_share(&dealing.commitments, 1, own_share, g), "kept share is valid");

        for e in &dealing.entries {
            let k = e.recipient;
            // Figure 4, line 9: g^z == R * X̄  (additively: z*g == R + X̄)
            let x_bar = share_commitment(&dealing.commitments, k);
            assert_eq!(g * e.z, e.big_r + x_bar, "public ciphertext check for party {k}");

            // Recipient k recomputes the pad cheaply and decrypts.
            let r = compute_r(&params, sks[(k - 1) as usize], &pks[(k - 1) as usize], &pks[0], &dealing.msg);
            assert_eq!(g * r, e.big_r, "compute_r must match the dealer's R");
            let share = e.z - r;
            assert!(verify_share(&dealing.commitments, k, share, g), "decrypted share must verify");
        }

        // The eVRF proof for one entry verifies (the other costs the same).
        let e = &dealing.entries[0];
        assert!(verify(&params, &pks[0], &pks[(e.recipient - 1) as usize], &dealing.msg, e.big_r, &e.proof));
    }

    #[test]
    fn compute_r_is_symmetric() {
        // Cheap: no proofs. Alice's and Bob's pads must agree.
        let params = setup_evrf();
        let (sks, pks) = parties(2);
        let r_ab = compute_r(&params, sks[0], &pks[0], &pks[1], b"m");
        let r_ba = compute_r(&params, sks[1], &pks[1], &pks[0], b"m");
        assert_eq!(r_ab, r_ba);
        assert_ne!(r_ab, compute_r(&params, sks[0], &pks[0], &pks[1], b"other"), "msg must matter");
    }
}
