#[cfg(test)]
mod tests {
    use crate::feldman::{deal, share_commitment, verify_share};

    use ark_bls12_381::{Fr, G1Projective};
    use ark_ec::PrimeGroup;
    use ark_ff::{Field, One, UniformRand, Zero};

    /// Lagrange interpolation at X = 0 from points (k, f(k)).
    fn reconstruct(points: &[(u64, Fr)]) -> Fr {
        let mut secret = Fr::zero();
        for (i, &(xi, yi)) in points.iter().enumerate() {
            let xi_f = Fr::from(xi);
            let mut lambda = Fr::one();
            for (j, &(xj, _)) in points.iter().enumerate() {
                if i != j {
                    let xj_f = Fr::from(xj);
                    lambda *= xj_f * (xj_f - xi_f).inverse().unwrap();
                }
            }
            secret += lambda * yi;
        }
        secret
    }

    fn setup(t: usize, n: usize) -> (Fr, G1Projective, Vec<G1Projective>, Vec<Fr>) {
        let mut rng = ark_std::test_rng();
        let g = G1Projective::generator();
        let secret = Fr::rand(&mut rng);
        let (commitments, shares) = deal(secret, t, n, g, &mut rng);
        (secret, g, commitments, shares)
    }

    #[test]
    fn sizes_match_threshold() {
        let (_, _, commitments, shares) = setup(3, 5);
        assert_eq!(commitments.len(), 3, "t commitments (degree t-1)");
        assert_eq!(shares.len(), 5, "one share per party");
    }

    #[test]
    fn honest_shares_verify() {
        let (_, g, commitments, shares) = setup(3, 5);
        for (i, s) in shares.iter().enumerate() {
            assert!(verify_share(&commitments, (i + 1) as u64, *s, g), "share {} must verify", i + 1);
        }
    }

    #[test]
    fn tampered_share_fails() {
        let (_, g, commitments, shares) = setup(3, 5);
        assert!(!verify_share(&commitments, 2, shares[1] + Fr::one(), g));
        // a correct share checked at the wrong index also fails
        assert!(!verify_share(&commitments, 3, shares[1], g));
    }

    #[test]
    fn a0_commits_to_the_secret() {
        let (secret, g, commitments, _) = setup(3, 5);
        assert_eq!(commitments[0], g * secret, "A0 = g^secret");
        assert_eq!(share_commitment(&commitments, 0), g * secret, "f(0) in the exponent");
    }

    #[test]
    fn share_commitment_is_g_to_the_share() {
        let (_, g, commitments, shares) = setup(4, 7);
        for (i, s) in shares.iter().enumerate() {
            assert_eq!(share_commitment(&commitments, (i + 1) as u64), g * s);
        }
    }

    #[test]
    fn any_t_shares_reconstruct_the_secret() {
        let (secret, _, _, shares) = setup(3, 5);
        let pts = |idx: &[u64]| -> Vec<(u64, Fr)> {
            idx.iter().map(|&k| (k, shares[(k - 1) as usize])).collect()
        };
        assert_eq!(reconstruct(&pts(&[1, 2, 3])), secret);
        assert_eq!(reconstruct(&pts(&[2, 4, 5])), secret);
        assert_eq!(reconstruct(&pts(&[1, 3, 5])), secret);
    }

    #[test]
    fn fewer_than_t_shares_do_not_reconstruct() {
        let (secret, _, _, shares) = setup(3, 5);
        let two = vec![(1u64, shares[0]), (2u64, shares[1])];
        assert_ne!(reconstruct(&two), secret, "t-1 shares must not give the secret");
    }
}
