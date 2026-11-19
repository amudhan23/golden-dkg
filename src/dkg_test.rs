#[cfg(test)]
mod tests {
    use crate::dkg::round0;
    use crate::evrf::{compute_r, setup_evrf, verify, EvrfParams};
    use crate::feldman::{share_commitment, verify_share};

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

    /// Lagrange interpolation at X = 0 from points (k, y_k).
    fn interpolate_at_zero(points: &[(u64, ark_bls12_381::Fr)]) -> ark_bls12_381::Fr {
        use ark_bls12_381::Fr;
        use ark_ff::{Field, One, Zero};
        let mut acc = Fr::zero();
        for (a, &(xa, ya)) in points.iter().enumerate() {
            let mut lambda = Fr::one();
            for (b, &(xb, _)) in points.iter().enumerate() {
                if a != b {
                    lambda *= Fr::from(xb) * (Fr::from(xb) - Fr::from(xa)).inverse().unwrap();
                }
            }
            acc += lambda * ya;
        }
        acc
    }

    /// The whole Golden DKG (Figure 4) with 3 parties and threshold 2.
    /// Slow (6 proofs + 12 verifications, ~3 min):
    ///     cargo test --release golden_dkg_end_to_end -- --ignored --nocapture
    #[test]
    #[ignore]
    fn golden_dkg_end_to_end_three_parties() {
        use crate::dkg::round1;
        use std::time::Instant;

        let params = setup_evrf();
        let (n, t) = (3usize, 2usize);
        let (sks, pks) = parties(n);
        let mut rng = ark_std::test_rng();
        let g = params.bp.u;

        // ── Round 0: every party deals ──
        let t0 = Instant::now();
        let mut dealings = Vec::new();
        let mut states = Vec::new();
        for i in 1..=n as u64 {
            let (d, st) = round0(&params, n, t, i, sks[(i - 1) as usize], &pks, &mut rng);
            dealings.push(d);
            states.push(st);
        }
        println!("round 0 (all parties): {:.1} s", t0.elapsed().as_secs_f64());

        // ── Round 1: every party verifies, decrypts, and outputs ──
        let t1 = Instant::now();
        let outputs: Vec<_> = (1..=n as u64)
            .map(|i| {
                round1(&params, n, i, sks[(i - 1) as usize], &pks, &dealings, states[(i - 1) as usize])
                    .expect("honest run must not abort")
            })
            .collect();
        println!("round 1 (all parties): {:.1} s", t1.elapsed().as_secs_f64());

        // Everyone agrees on the group key and on every party's public key share.
        for o in &outputs[1..] {
            assert_eq!(o.pk, outputs[0].pk, "all parties must output the same PK");
            assert_eq!(o.pk_shares, outputs[0].pk_shares, "all parties must agree on PK_l");
        }

        // Each party's secret share matches its public key share: g * sk_l == PK_l.
        for (l, o) in outputs.iter().enumerate() {
            assert_eq!(g * o.sk, outputs[0].pk_shares[l], "sk_{} must match PK_{}", l + 1, l + 1);
        }

        // Any t parties reconstruct the secret behind PK; nobody ever held it alone.
        for pair in [[1u64, 2], [1, 3], [2, 3]] {
            let pts: Vec<_> = pair.iter().map(|&k| (k, outputs[(k - 1) as usize].sk)).collect();
            let secret = interpolate_at_zero(&pts);
            assert_eq!(g * secret, outputs[0].pk, "parties {:?} must reconstruct the key behind PK", pair);
        }
        // t-1 = 1 party alone does not.
        assert_ne!(g * outputs[0].sk, outputs[0].pk, "one share alone is not the group secret");

        // ── A cheating dealer is caught ──
        // Dealer 1 corrupts the ciphertext it sent to party 2.
        dealings[0].entries[0].z += ark_bls12_381::Fr::from(1u64);
        assert!(
            round1(&params, n, 3, sks[2], &pks, &dealings, states[2]).is_none(),
            "a bad ciphertext must make honest parties abort"
        );
    }
}
