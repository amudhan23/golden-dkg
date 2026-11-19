use ark_bls12_381::G1Projective;
use ark_ec::AffineRepr;
use ark_ec::CurveGroup;
use ark_ec::twisted_edwards::TECurveConfig;
use ark_ed_on_bls12_381_bandersnatch::BandersnatchConfig;
use ark_ed_on_bls12_381_bandersnatch::EdwardsAffine;
use ark_ed_on_bls12_381_bandersnatch::EdwardsProjective;
use ark_ed_on_bls12_381_bandersnatch::Fq;
use ark_ed_on_bls12_381_bandersnatch::Fr as Scalar;
use ark_ff::PrimeField;
use ark_ff::{One, Zero};
use ark_relations::r1cs::ConstraintSystem;
use ark_relations::r1cs::SynthesisMode;
use ark_relations::r1cs::Variable;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::Rng;
use bulletproofs_ipa::ProofParams;
use bulletproofs_ipa::protocol;
use merlin::Transcript;
use sha2::{Digest, Sha256};

use crate::evrf_dh_circuit::evrf_dh;
use crate::schnorr;
use crate::schnorr::SchnorrProof; // your trait: append_point, append_scalar, get_challenge

fn hash_to_curve<G: CurveGroup>(domain: &[u8], data: &[u8]) -> G {
    for counter in 0u64.. {
        let mut hasher = Sha256::new();
        hasher.update(domain);
        hasher.update(b"/");
        hasher.update(data);
        hasher.update(b"/");
        hasher.update(counter.to_le_bytes());
        let hash = hasher.finalize();

        // Try to map these bytes to a curve point
        if let Some(point) = G::Affine::from_random_bytes(&hash) {
            let point = point.mul_by_cofactor_to_group();
            if !point.is_zero() {
                return point;
            }
        }
    }
    unreachable!()
}

pub fn hash_points(
    msg: &[u8],
    pk_a: &EdwardsAffine,
    pk_b: &EdwardsAffine,
) -> (EdwardsAffine, EdwardsAffine) {
    let mut pk_a_vec = Vec::new();
    pk_a.serialize_compressed(&mut pk_a_vec).unwrap();

    let mut pk_b_vec = Vec::new();
    pk_b.serialize_compressed(&mut pk_b_vec).unwrap();

    let (first, second) = if pk_a_vec <= pk_b_vec {
        (pk_a_vec, pk_b_vec)
    } else {
        (pk_b_vec, pk_a_vec)
    };

    let mut data = msg.to_vec();
    data.extend_from_slice(&first);
    data.extend_from_slice(&second);

    let h1 = hash_to_curve::<EdwardsProjective>(b"golden_dkg/evrf/H1", &data).into_affine();
    let h2 = hash_to_curve::<EdwardsProjective>(b"golden_dkg/evrf/H2", &data).into_affine();

    (h1, h2)
}

pub struct EvrfParams {
    pub bp: ProofParams<G1Projective>,
    pub g_in: EdwardsAffine,
    pub a: Fq,
    pub d: Fq,
    pub beta: Fq,
}

pub fn setup_evrf() -> EvrfParams {
    let cs = ConstraintSystem::<Fq>::new_ref();
    cs.set_mode(SynthesisMode::Setup);
    evrf_dh(
        cs.clone(),
        Fq::zero(),
        Fq::one(), // g_in      (dummy point: the identity (0, 1))
        Fq::zero(),
        Fq::one(), // PK₁
        Fq::zero(),
        Fq::one(), // PK₂
        Fq::zero(),
        Fq::one(), // H₁
        Fq::zero(),
        Fq::one(),  // H₂
        Fq::zero(), // β
        Fq::zero(), // sk (dummy)
        BandersnatchConfig::COEFF_A,
        BandersnatchConfig::COEFF_D,
    )
    .unwrap();
    let n = (cs.num_constraints() + cs.num_witness_variables()).next_power_of_two();
    let bp = bulletproofs_ipa::setup_params(n);
    let g_in = EdwardsAffine::generator();
    let a = BandersnatchConfig::COEFF_A;
    let d = BandersnatchConfig::COEFF_D;
    let beta = Fq::from_le_bytes_mod_order(&Sha256::digest(b"golden_dkg/evrf/beta"));

    EvrfParams {
        bp,
        g_in,
        a,
        d,
        beta,
    }
}

pub struct EvrfProof {
    pub bp: bulletproofs_ipa::protocol::Proof<G1Projective>,
    pub pok: SchnorrProof<G1Projective>,
}

pub fn evaluate(
    params: &EvrfParams,
    sk: Scalar,
    my_pk: &EdwardsAffine,
    other_pk: &EdwardsAffine,
    msg: &[u8],
    rng: &mut impl Rng,
) -> (Fq, G1Projective, EvrfProof) {
    let (h1, h2) = hash_points(msg, my_pk, other_pk);
    let sk_fq = Fq::from_bigint(sk.into_bigint().into()).unwrap();
    let cs = ConstraintSystem::<Fq>::new_ref();

    let (r_var, r) = evrf_dh(
        cs.clone(),
        params.g_in.x,
        params.g_in.y,
        my_pk.x,
        my_pk.y,
        other_pk.x,
        other_pk.y,
        h1.x,
        h1.y,
        h2.x,
        h2.y,
        params.beta,
        sk_fq,
        params.a,
        params.d,
    )
    .unwrap();

    cs.finalize();

    let mats = cs.to_matrices().unwrap();
    let b = cs.borrow().unwrap();
    let instance = b.instance_assignment.clone();
    let witness = b.witness_assignment.clone();
    drop(b);

    let r_index = match r_var {
        Variable::Witness(i) => i,
        _ => panic!("r must be a witness variable"),
    };

    let mut transcript = Transcript::new(b"unit-golden");

    let bp_proof = bulletproofs_ipa::protocol::prove(
        &mats,
        &instance,
        &witness,
        &params.bp,
        rng,
        &[r_index],
        &[Fq::zero()],
        &mut transcript,
    );

    let R = params.bp.u * r;

    let pok = schnorr::prove(params.bp.u, R, r, rng, &mut transcript);

    (r, R, EvrfProof { bp: bp_proof, pok })
}

// pub fn evrf_dh<F: PrimeField>(
//     cs: ConstraintSystemRef<F>,
//     g_in_x: F,
//     g_in_y: F,
//     pk1_x: F,
//     pk1_y: F,
//     pk2_x: F,
//     pk2_y: F,
//     h1_x: F,
//     h1_y: F,
//     h2_x: F,
//     h2_y: F,
//     beta: F,
//     witness: F,
//     a: F,
//     d: F,
// )
