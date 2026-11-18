use ark_ec::CurveGroup;
use ark_ff::UniformRand;
use ark_std::rand::Rng;
use bulletproofs_ipa::TranscriptProtocol; // your trait: append_point, append_scalar, get_challenge
use merlin::Transcript;

pub struct SchnorrProof<G: CurveGroup> {
    pub A: G,
    pub s: G::ScalarField,
}

pub fn prove<G: CurveGroup>(
    base: G,
    public: G,
    secret: G::ScalarField,
    rng: &mut impl Rng,
    transcript: &mut Transcript,
) -> SchnorrProof<G> {
    let k = G::ScalarField::rand(rng);
    let A = base * k;
    transcript.append_point(b"base", &base);
    transcript.append_point(b"public", &public);
    transcript.append_point(b"A", &A);
    let c: G::ScalarField = transcript.get_challenge(b"c");
    let s = k + c * secret;

    SchnorrProof { A, s }
}

pub fn verify<G: CurveGroup>(
    schnorr_proof: &SchnorrProof<G>,
    base: G,
    public: G,
    transcript: &mut Transcript,
) -> bool {
    transcript.append_point(b"base", &base);
    transcript.append_point(b"public", &public);
    transcript.append_point(b"A", &schnorr_proof.A);
    let c: G::ScalarField = transcript.get_challenge(b"c");
    base * schnorr_proof.s == schnorr_proof.A + public * c
}
