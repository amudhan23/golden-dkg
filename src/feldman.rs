use ark_bls12_381::{Fr, G1Projective};
use ark_ec::{AffineRepr, CurveGroup};
use ark_ed_on_bls12_381_bandersnatch::{EdwardsAffine, Fq, Fr as Scalar};
use ark_ff::{One, PrimeField, UniformRand, Zero};
use ark_std::rand::Rng;

use crate::evrf::{EvrfParams, hash_points};

pub fn deal(
    secret: Fr,
    t: usize,
    n: usize,
    g: G1Projective,
    rng: &mut impl Rng,
) -> (Vec<G1Projective>, Vec<Fr>) {
    let mut poly_coeff = Vec::new();
    let mut shares = Vec::new();

    poly_coeff.push(secret);

    for _ in 1..t {
        let coeff = Fr::rand(rng);
        poly_coeff.push(coeff);
    }

    let commitments: Vec<G1Projective> = poly_coeff.iter().map(|a| g * a).collect();

    for i in 1..n + 1 {
        let mut share = Fr::zero();

        let x = Fr::from(i as u64);
        for l in (0..poly_coeff.len()).rev() {
            share = share * x + poly_coeff[l];
        }
        shares.push(share);
    }

    (commitments, shares)
}

pub fn share_commitment(commitments: &[G1Projective], k: u64) -> G1Projective {
    let x = Fr::from(k as u64);

    let mut acc = G1Projective::zero();
    for i in (0..commitments.len()).rev() {
        acc = acc * x + commitments[i];
    }

    acc
}

pub fn verify_share(commitments: &[G1Projective], k: u64, share: Fr, g: G1Projective) -> bool {
    g * share == share_commitment(commitments, k)
}
