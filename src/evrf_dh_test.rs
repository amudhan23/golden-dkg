#[cfg(test)]
mod tests {
    use crate::evrf_dh_circuit::evrf_dh;

    use ark_ec::twisted_edwards::TECurveConfig;
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ed_on_bls12_381_bandersnatch::{BandersnatchConfig, EdwardsAffine, Fq, Fr as Scalar};
    use ark_ff::{One, PrimeField};
    use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef, Variable};
    use ark_std::UniformRand;

    /// Everything the relation R_T-eVRF-DH needs, computed natively (outside the circuit).
    struct Instance {
        g: EdwardsAffine,
        pk1: EdwardsAffine,
        pk2: EdwardsAffine,
        h1: EdwardsAffine,
        h2: EdwardsAffine,
        beta: Fq,
        sk1_fq: Fq,
        expected_r: Fq,
    }

    /// Embed a Bandersnatch scalar (F_s) into the base field (F_p). Works because s < p.
    fn scalar_to_fq(s: Scalar) -> Fq {
        Fq::from_bigint(s.into_bigint().into()).unwrap()
    }

    fn random_point<R: ark_std::rand::Rng>(rng: &mut R) -> EdwardsAffine {
        (EdwardsAffine::generator() * Scalar::rand(rng)).into_affine()
    }

    fn make_instance() -> Instance {
        let mut rng = ark_std::test_rng();
        let g = EdwardsAffine::generator();

        // (0) Party 1's key pair
        let sk1 = Scalar::rand(&mut rng);
        let pk1 = (g * sk1).into_affine();

        // Party 2's public key (party 2's secret is never used by party 1)
        let pk2 = random_point(&mut rng);

        // Hash-to-curve outputs H1, H2. For this test, any two independent
        // subgroup points stand in for H_{G_in,j}(msg, PK1, PK2).
        let h1 = random_point(&mut rng);
        let h2 = random_point(&mut rng);
        let beta = Fq::rand(&mut rng);

        // (1) S = PK2^sk1
        let s = (pk2 * sk1).into_affine();
        // (2)-(3) k = int(S.x): S.x is an F_p element, used as an integer exponent
        let k = s.x.into_bigint();
        // (4)-(5) T_j = H_j^k (integer scalar multiplication, same as the circuit's bits)
        let t1 = h1.mul_bigint(k).into_affine();
        let t2 = h2.mul_bigint(k).into_affine();
        // (6)-(8) r = beta * T1.x + T2.x
        let expected_r = beta * t1.x + t2.x;

        Instance {
            g,
            pk1,
            pk2,
            h1,
            h2,
            beta,
            sk1_fq: scalar_to_fq(sk1),
            expected_r,
        }
    }

    /// Run the circuit on `inst`, optionally overriding the claimed PK1 or the secret key.
    fn run(
        inst: &Instance,
        pk1: EdwardsAffine,
        sk1_fq: Fq,
    ) -> (ConstraintSystemRef<Fq>, Variable, Fq) {
        let cs = ConstraintSystem::<Fq>::new_ref();
        let a: Fq = BandersnatchConfig::COEFF_A;
        let d: Fq = BandersnatchConfig::COEFF_D;

        let (r_var, r) = evrf_dh(
            cs.clone(),
            inst.g.x,
            inst.g.y,
            pk1.x,
            pk1.y,
            inst.pk2.x,
            inst.pk2.y,
            inst.h1.x,
            inst.h1.y,
            inst.h2.x,
            inst.h2.y,
            inst.beta,
            sk1_fq,
            a,
            d,
        )
        .unwrap();
        (cs, r_var, r)
    }

    #[test]
    fn evrf_dh_honest_matches_native() {
        let inst = make_instance();
        let (cs, _, r) = run(&inst, inst.pk1, inst.sk1_fq);

        assert_eq!(r, inst.expected_r, "circuit r doesn't match native r");
        assert!(cs.is_satisfied().unwrap(), "honest eVRF should satisfy all constraints");

        println!(
            "evrf_dh: {} constraints, {} witness vars, {} public inputs",
            cs.num_constraints(),
            cs.num_witness_variables(),
            cs.num_instance_variables()
        );
    }

    #[test]
    fn evrf_dh_wrong_pk1_fails() {
        // Claimed PK1 isn't g^sk1 -> constraint (0) must fail.
        let inst = make_instance();
        let wrong_pk1 = (inst.pk1 + inst.g).into_affine();
        assert_ne!(wrong_pk1, inst.pk1);

        let (cs, _, _) = run(&inst, wrong_pk1, inst.sk1_fq);
        assert!(!cs.is_satisfied().unwrap(), "wrong PK1 must not satisfy");
    }

    #[test]
    fn evrf_dh_wrong_sk_fails() {
        // Prover doesn't know the secret key for PK1 -> constraint (0) must fail.
        let inst = make_instance();
        let wrong_sk = inst.sk1_fq + Fq::one();

        let (cs, _, _) = run(&inst, inst.pk1, wrong_sk);
        assert!(!cs.is_satisfied().unwrap(), "wrong sk1 must not satisfy");
    }

    #[test]
    fn evrf_dh_tampered_r_fails() {
        // Honest run, then overwrite r's witness value -> constraint (8) must fail.
        let inst = make_instance();
        let (cs, r_var, _) = run(&inst, inst.pk1, inst.sk1_fq);
        assert!(cs.is_satisfied().unwrap());

        match r_var {
            Variable::Witness(i) => {
                cs.borrow_mut().unwrap().witness_assignment[i] += Fq::one();
            }
            _ => panic!("r should be a witness variable"),
        }
        assert!(!cs.is_satisfied().unwrap(), "tampered r must not satisfy");
    }
}
