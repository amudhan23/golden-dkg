#[cfg(test)]
mod tests {
    use crate::schnorr::{prove, verify, SchnorrProof};

    use ark_ec::{CurveGroup, PrimeGroup};
    use ark_ff::{Field, One, UniformRand};
    use bulletproofs_ipa::TranscriptProtocol;
    use merlin::Transcript;

    const LABEL: &[u8] = b"golden_dkg/schnorr_test";

    /// Honest prove + verify for any curve.
    fn roundtrip<G: CurveGroup + PrimeGroup>() {
        let mut rng = ark_std::test_rng();
        let base = G::generator();
        let x = G::ScalarField::rand(&mut rng);
        let public = base * x;

        let proof = prove(base, public, x, &mut rng, &mut Transcript::new(LABEL));
        assert!(verify(&proof, base, public, &mut Transcript::new(LABEL)));
    }

    #[test]
    fn honest_verifies_on_bandersnatch() {
        // identity keys: PK = sk * g_in
        roundtrip::<ark_ed_on_bls12_381_bandersnatch::EdwardsProjective>();
    }

    #[test]
    fn honest_verifies_on_bls12_381() {
        // eVRF output commitment: R = r * g_out
        roundtrip::<ark_bls12_381::G1Projective>();
    }

    type G = ark_ed_on_bls12_381_bandersnatch::EdwardsProjective;
    type F = <G as PrimeGroup>::ScalarField;

    fn setup() -> (G, F, G, SchnorrProof<G>) {
        let mut rng = ark_std::test_rng();
        let base = G::generator();
        let x = F::rand(&mut rng);
        let public = base * x;
        let proof = prove(base, public, x, &mut rng, &mut Transcript::new(LABEL));
        (base, x, public, proof)
    }

    #[test]
    fn wrong_public_key_fails() {
        let (base, _, public, proof) = setup();
        assert!(!verify(&proof, base, public + base, &mut Transcript::new(LABEL)));
    }

    #[test]
    fn wrong_base_fails() {
        let (base, _, public, proof) = setup();
        assert!(!verify(&proof, base + base, public, &mut Transcript::new(LABEL)));
    }

    #[test]
    fn tampered_s_fails() {
        let (base, _, public, mut proof) = setup();
        proof.s += F::one();
        assert!(!verify(&proof, base, public, &mut Transcript::new(LABEL)));
    }

    #[test]
    fn tampered_a_fails() {
        let (base, _, public, mut proof) = setup();
        proof.A += base;
        assert!(!verify(&proof, base, public, &mut Transcript::new(LABEL)));
    }

    #[test]
    fn different_transcript_label_fails() {
        let (base, _, public, proof) = setup();
        assert!(!verify(&proof, base, public, &mut Transcript::new(b"some other context")));
    }

    #[test]
    fn prover_without_secret_fails() {
        // A cheater who doesn't know x uses a guess instead.
        let mut rng = ark_std::test_rng();
        let base = G::generator();
        let public = base * F::rand(&mut rng);
        let guess = F::rand(&mut rng);
        let proof = prove(base, public, guess, &mut rng, &mut Transcript::new(LABEL));
        assert!(!verify(&proof, base, public, &mut Transcript::new(LABEL)));
    }

    #[test]
    fn statement_is_bound_into_challenge() {
        // Frozen Heart check: choosing `public` AFTER seeing the challenge must fail.
        // Forger picks A and s freely, computes c the way a WEAK transcript would
        // (from A alone), solves for a public key X = (s*B - A)/c, and submits.
        let mut rng = ark_std::test_rng();
        let base = G::generator();
        let a = base * F::rand(&mut rng);
        let s = F::rand(&mut rng);
        let mut weak = Transcript::new(LABEL);
        weak.append_point(b"A", &a);
        let c: F = weak.get_challenge(b"c");
        let forged_public = (base * s - a) * c.inverse().unwrap();
        // It satisfies the equation for the weak challenge...
        assert_eq!(base * s, a + forged_public * c);
        // ...but our verifier hashes base and public too, so the challenge differs.
        let proof = SchnorrProof { A: a, s };
        assert!(!verify(&proof, base, forged_public, &mut Transcript::new(LABEL)));
    }

    #[test]
    fn nonce_reuse_leaks_the_secret() {
        // Special soundness: two responses for the same A (same k) give x.
        let mut rng = ark_std::test_rng();
        let base = G::generator();
        let x = F::rand(&mut rng);
        let k = F::rand(&mut rng);
        let (c1, c2) = (F::rand(&mut rng), F::rand(&mut rng));
        let s1 = k + c1 * x;
        let s2 = k + c2 * x;
        let recovered = (s1 - s2) * (c1 - c2).inverse().unwrap();
        assert_eq!(recovered, x, "x = (s1 - s2) / (c1 - c2)");
        let _ = base;
    }
}
