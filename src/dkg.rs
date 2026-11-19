use ark_bls12_381::{Fr, G1Projective};
use ark_ed_on_bls12_381_bandersnatch::EdwardsAffine;
use ark_ed_on_bls12_381_bandersnatch::Fr as Scalar;
use ark_ff::UniformRand as _;
use ark_std::rand;
use rand::Rng;

use crate::evrf::evaluate;
use crate::evrf::{EvrfParams, EvrfProof};
use crate::feldman;

pub struct RecipientEntry {
    pub recipient: u64,
    pub big_r: G1Projective,
    pub z: Fr,
    pub proof: EvrfProof,
}

pub struct Dealing {
    // bmsg_i
    pub dealer: u64,
    pub msg: Vec<u8>,
    pub commitments: Vec<G1Projective>, // C̄_i = (A_0 … A_{t-1})
    pub entries: Vec<RecipientEntry>,   // one per j ≠ i
}

pub struct DkgOutput {
    pub pk: G1Projective,
    pub pk_shares: Vec<G1Projective>,
    pub sk: Fr,
}

pub fn round0(
    params: &EvrfParams,
    n: usize,
    t: usize,
    i: u64,
    sk_i: Scalar,
    pks: &[EdwardsAffine],
    rng: &mut impl Rng,
) -> (Dealing, Fr) {
    let omega = Fr::rand(rng);
    let (commitments, shares) = feldman::deal(omega, t, n, params.bp.u, rng);
    let msg: Vec<u8> = (0..32).map(|_| rng.r#gen::<u8>()).collect();
    let mut recipientEntries = Vec::new();

    let i_us = i as usize;
    for j in 1..n + 1 {
        if j == i_us {
            continue;
        }

        let (r, big_r, proof) = evaluate(params, sk_i, &pks[i_us - 1], &pks[j - 1], &msg, rng);
        let z = r + shares[j - 1];
        let recipientEntry = RecipientEntry {
            recipient: j as u64,
            big_r,
            z,
            proof,
        };
        recipientEntries.push(recipientEntry);
    }

    (
        Dealing {
            dealer: i,
            msg,
            commitments,
            entries: recipientEntries,
        },
        shares[i_us - 1],
    )
}
