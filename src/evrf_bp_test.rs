//! Prove and verify the eVRF-DH circuit with our own Bulletproofs.
//!
//! The prover and the verifier each run `evrf_dh` themselves:
//! - prover:   Prove mode, real sk1  -> matrices + instance + witness
//! - verifier: Setup mode, dummy sk1 -> its OWN matrices (never the prover's)
//!
//! Slow (n = 65,536 gates), so it is `#[ignore]`d. Run with:
//!     cargo test --release evrf_bulletproof -- --ignored --nocapture

#[cfg(test)]
mod tests {
    use crate::evrf_dh_circuit::evrf_dh;

    use ark_bls12_381::G1Projective;
    use ark_ec::twisted_edwards::TECurveConfig;
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ed_on_bls12_381_bandersnatch::{BandersnatchConfig, EdwardsAffine, Fq, Fr as Scalar};
    use ark_ff::{One, PrimeField, Zero};
    use ark_relations::r1cs::{ConstraintMatrices, ConstraintSystem, SynthesisMode, Variable};
    use ark_std::UniformRand;
    use bulletproofs_ipa::protocol::{prove, verify};
    use bulletproofs_ipa::setup_params;
    use merlin::Transcript;
    use std::time::Instant;

    const LABEL: &[u8] = b"golden_dkg/evrf_dh/v1";

    /// Public statement of R_T-eVRF-DH (everything the verifier knows).
    struct Public {
        g: EdwardsAffine,
        pk1: EdwardsAffine,
        pk2: EdwardsAffine,
        h1: EdwardsAffine,
        h2: EdwardsAffine,
        beta: Fq,
    }

    impl Public {
        /// The instance vector in the exact order `evrf_dh` allocates its
        /// public inputs: [1, g, PK1, PK2, H1, H2, beta].
        fn instance(&self) -> Vec<Fq> {
            vec![
                Fq::one(),
                self.g.x, self.g.y,
                self.pk1.x, self.pk1.y,
                self.pk2.x, self.pk2.y,
                self.h1.x, self.h1.y,
                self.h2.x, self.h2.y,
                self.beta,
            ]
        }
    }

    /// Run the circuit in the given mode and return the matrices (+ assignments).
    fn synthesize(
        public: &Public,
        sk1: Fq,
        mode: SynthesisMode,
    ) -> (ConstraintMatrices<Fq>, Vec<Fq>, Vec<Fq>, usize) {
        let cs = ConstraintSystem::<Fq>::new_ref();
        cs.set_mode(mode);
        let a: Fq = BandersnatchConfig::COEFF_A;
        let d: Fq = BandersnatchConfig::COEFF_D;
        let (r_var, _) = evrf_dh(
            cs.clone(),
            public.g.x, public.g.y,
            public.pk1.x, public.pk1.y,
            public.pk2.x, public.pk2.y,
            public.h1.x, public.h1.y,
            public.h2.x, public.h2.y,
            public.beta,
            sk1,
            a, d,
        )
        .unwrap();
        cs.finalize();
        let mats = cs.to_matrices().unwrap();
        let b = cs.borrow().unwrap();
        let r_index = match r_var {
            Variable::Witness(i) => i,
            _ => panic!("r must be a witness variable"),
        };
        (mats, b.instance_assignment.clone(), b.witness_assignment.clone(), r_index)
    }

    #[test]
    #[ignore]
    fn evrf_bulletproof_prove_and_verify() {
        let mut rng = ark_std::test_rng();
        let g = EdwardsAffine::generator();

        // ── Statement and witness ──────────────────────────────────────
        let sk1 = Scalar::rand(&mut rng);
        let sk1_fq = Fq::from_bigint(sk1.into_bigint().into()).unwrap();
        let rand_point = |rng: &mut _| (g * Scalar::rand(rng)).into_affine();
        let public = Public {
            g,
            pk1: (g * sk1).into_affine(),
            pk2: rand_point(&mut rng),
            h1: rand_point(&mut rng),
            h2: rand_point(&mut rng),
            beta: Fq::rand(&mut rng),
        };

        // ── Prover: Prove mode, real secret ────────────────────────────
        let (mats_p, instance_p, witness, r_index_p) =
            synthesize(&public, sk1_fq, SynthesisMode::Prove { construct_matrices: true });

        // ── Verifier: Setup mode, dummy secret, its own matrices ───────
        let (mats_v, instance_setup, witness_setup, r_index_v) =
            synthesize(&public, Fq::zero(), SynthesisMode::Setup);
        let instance_v = public.instance();

        assert!(witness_setup.is_empty(), "Setup mode must not record witness values");
        assert_eq!(mats_v.a, mats_p.a, "verifier's A must equal prover's A");
        assert_eq!(mats_v.b, mats_p.b, "verifier's B must equal prover's B");
        assert_eq!(mats_v.c, mats_p.c, "verifier's C must equal prover's C");
        assert_eq!(instance_v, instance_p, "verifier's instance must equal prover's");
        assert_eq!(r_index_v, r_index_p, "verifier finds r at the same witness index");
        println!("setup-mode instance_assignment length: {}", instance_setup.len());

        let n = (mats_v.num_constraints + mats_v.num_witness_variables).next_power_of_two();
        println!(
            "constraints: {}, witness vars: {}, public inputs: {}, gates (padded): {}",
            mats_v.num_constraints,
            mats_v.num_witness_variables,
            mats_v.num_instance_variables,
            n
        );

        // ── Public parameters (deterministic, both sides derive them) ──
        let t = Instant::now();
        let params = setup_params::<G1Projective>(n);
        assert_eq!(params.g.len(), n);
        assert_eq!(params.h.len(), n);
        println!("setup_params:  {:.2} s", t.elapsed().as_secs_f64());

        // ── Prove ──────────────────────────────────────────────────────
        let t = Instant::now();
        let mut tp = Transcript::new(LABEL);
        let r = witness[r_index_p];
        let proof = prove::<G1Projective>(
            &mats_p, &instance_p, &witness, &params, &mut rng,
            &[r_index_p], &[Fq::zero()], &mut tp,
        );
        // What the dealer posts on the board: R = r * g_out, with g_out = U.
        let big_r = params.u * r;
        println!("prove:         {:.2} s", t.elapsed().as_secs_f64());

        // ── Verify (with the verifier's own matrices and instance) ─────
        let t = Instant::now();
        let mut tv = Transcript::new(LABEL);
        let ok = verify::<G1Projective>(
            &mats_v, &instance_v, &params, &[r_index_v], &[big_r], &proof, &mut tv,
        );
        println!("verify:        {:.2} s", t.elapsed().as_secs_f64());
        assert!(ok, "honest eVRF proof must verify");

        // ── Proof size: 8 points + 3 scalars + IPA (2 log n points + 2 scalars)
        let rounds = n.trailing_zeros() as usize;
        let points = 8 + 2 * rounds; // R is posted separately, not part of the proof
        let scalars = 3 + 2;
        println!(
            "proof size:    {} G1 points + {} scalars = {} bytes (compressed)",
            points,
            scalars,
            points * 48 + scalars * 32
        );

        // ── Negative: same proof, verifier claims a different PK1 ──────
        let mut wrong = public.instance();
        let wrong_pk1 = (public.pk1 + public.g).into_affine();
        wrong[3] = wrong_pk1.x;
        wrong[4] = wrong_pk1.y;
        let mut tv = Transcript::new(LABEL);
        assert!(
            !verify::<G1Projective>(&mats_v, &wrong, &params, &[r_index_v], &[big_r], &proof, &mut tv),
            "proof must not verify against a different PK1"
        );

        // ── Negative: dealer posts R for a different r (the fake-R attack)
        let fake_r = params.u * (r + Fq::one());
        let mut tv = Transcript::new(LABEL);
        assert!(
            !verify::<G1Projective>(&mats_v, &instance_v, &params, &[r_index_v], &[fake_r], &proof, &mut tv),
            "proof must not verify against an R for a different r"
        );

        // ── Negative: same proof, different beta ───────────────────────
        let mut wrong = public.instance();
        wrong[11] += Fq::one();
        let mut tv = Transcript::new(LABEL);
        assert!(
            !verify::<G1Projective>(&mats_v, &wrong, &params, &[r_index_v], &[big_r], &proof, &mut tv),
            "proof must not verify against a different beta"
        );
    }
}
