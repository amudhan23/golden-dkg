# golden-dkg

A Rust implementation of **Golden**, the one-round distributed key generation
protocol of Bünz, Choi and Komlo ([ePrint 2025/1924](https://eprint.iacr.org/2025/1924)),
using the proof system the paper specifies: Bulletproofs over BLS12-381, with
the two-party eVRF proven as an R1CS over Bandersnatch.

The Bulletproofs prover/verifier is my own implementation
([bulletProofRust](https://github.com/amudhan23/bulletProofRust)) and is pulled
in as a git dependency.

## Status

Working end to end:

- the two-party eVRF (`Setup` / `Evaluate` / `Verify`, paper §4.2), with the
  output bound to `R = g_out^r` through a committed input, so step (9) of
  R_T-eVRF-DH costs no constraints (§5.4)
- Feldman VSS, share encryption `z = r + share`, and the public check
  `g^z = R · X̄`
- the full protocol of Figure 4 (Round 0 and Round 1)
- a 3-party, threshold-2 run in which all parties agree on the group key,
  every secret share matches its public share, any two parties reconstruct
  the secret, and a dealer that corrupts a ciphertext causes honest parties
  to abort

Not done yet: the batched variant of §6.2 (one proof per dealer instead of
n − 1), and the items under *Known limitations*.

## Results

eVRF proof for one statement, measured on a laptop with 16 threads, compared
with the paper's zkalc estimates (§5.5, Table 4, 8 vCPUs):

|                   | Paper (estimate) | This implementation |
|-------------------|------------------|---------------------|
| Constraints       | 3,598            | 3,585               |
| Vector length *n* | 8,192            | 8,192               |
| Prove             | 0.3 s            | 1.52 s              |
| Verify            | 0.1 s            | 0.24 s              |
| Proof size        | ~1.5 KB          | 1,792 B (+ 80 B Schnorr) |

The paper has no implementation; its numbers are estimates. The prove time
here is the whole `evaluate` call (circuit synthesis including the
precomputed tables, R1CS conversion, the Bulletproof, and the Schnorr proof),
while the estimate covers the dominant MSMs only. The proof is a little larger
because my Bulletproofs variant sends 34 group elements per proof against the
29 of the paper's Theorem 1.

How the numbers moved during development (2-thread machine):

| Version                                   | Constraints | *n*    | Prove  | Verify |
|-------------------------------------------|-------------|--------|--------|--------|
| Edwards double-and-add, naive MSM         | 16,909      | 65,536 | 120 s  | 69 s   |
| + Pippenger MSM, parallel                 | 16,909      | 65,536 | 54 s   | 50 s   |
| + precomputed public bases (Edwards)      | 5,645       | 16,384 | 10.2 s | 9.2 s  |
| + Weierstrass chord rule (§5.3)           | 3,585       | 8,192  | 5.6 s  | 5.5 s  |
| + single-MSM IPA verifier, parallel folding | 3,585     | 8,192  | 3.4 s  | 0.39 s |

A full 3-party DKG takes about 20 s for Round 0 and 5 s for Round 1 on the
same 2-thread machine, for all parties together.

## Design notes

**Two curves.** The DKG keys, Feldman commitments and the Bulletproofs
commitments live in BLS12-381 G1. Constraints are therefore over its scalar
field Fr. The eVRF has to compute Diffie–Hellman inside the circuit, which is
only cheap on a curve whose coordinates are in Fr, so identity keys live on
Bandersnatch. The eVRF output `r = β·T₁.x + T₂.x` is built from Bandersnatch
coordinates, which makes it an Fr element, so the same `r` can pad an Fr share
and be committed as `r·g_out`.

**Exponentiation gadget.** All four bases in R_T-eVRF-DH (g_in, PK₂, H₁, H₂)
are public, so their multiples 2ʲ·base are precomputed outside the circuit and
the bases end up as matrix coefficients rather than public inputs. Two gadgets
are implemented:

- `scalar_mul_fixed` (twisted Edwards): 5 constraints per bit, no correction
  points, because the selected point at bit 0 is exactly the identity and the
  addition law has no exceptional cases on the prime-order subgroup.
- `scalar_mul_sw_fixed` (short Weierstrass, paper §5.3): 3 constraints per bit
  with correction points. This is the one the eVRF uses.

**Correction points.** Two differences from §5.3:

- The last correction exponent is `c_λ = s − 2λ + 1`. With the paper's
  `s − 2λ − 1` the exponents sum to `s − 2`, so the corrections don't cancel.
- The corrections use an independent, hash-derived generator `h_c` instead of
  g. With g, step (0) multiplies g itself, the partial sums have small known
  exponents, and some steps could land on equal x-coordinates.

**Bandersnatch is not complete.** a = −5 is a non-square in Fr, so the
twisted Edwards addition law has exceptional pairs on the full curve. They
occur only when P ± Q is a point at infinity of order 2, which can't happen
for points in the odd-order subgroup. Verification therefore subgroup-checks
both public keys; H₁ and H₂ come from hash-to-curve with cofactor clearing.

**Committed inputs.** The Bulletproofs library takes a list of witness
indices and external Pedersen commitments `V = v·U + γ·H_b`. Each committed
value adds one linear row tying the wire to the commitment. With `γ = 0` and
`U = g_out`, the commitment is exactly the posted `R`.

## Layout

```
src/
  circuit_evrf.rs      bit decomposition
  add_gadget.rs        twisted Edwards addition gadgets
  scalar_mul.rs        Edwards scalar multiplication (double-and-add, fixed base)
  sw.rs                Weierstrass conversion, native arithmetic, §5.3 gadget
  evrf_dh_circuit.rs   R_T-eVRF-DH (Edwards and Weierstrass versions)
  evrf.rs              eVRF Setup / Evaluate / Verify, compute_r, hash to curve
  schnorr.rs           Schnorr proof of knowledge (generic over the curve)
  feldman.rs           Feldman VSS
  dkg.rs               Figure 4: round0, round1
  *_test.rs            tests for each module
```

## Running

```
cargo test --release                                   # fast tests (~5 s)

# slower tests, run explicitly:
cargo test --release evrf_sw_benchmark -- --ignored --nocapture
cargo test --release golden_dkg_end_to_end -- --ignored --nocapture
cargo test --release evaluate_both_sides -- --ignored --nocapture
```

## Known limitations

- **Fiat–Shamir does not bind the statement.** The Bulletproofs transcript
  absorbs the prover's messages and the committed inputs but not the public
  instance or the circuit. A prover able to pick a public input after seeing
  the challenges could satisfy the combined check for an unsatisfied circuit.
  Fix: absorb the instance and a hash of the matrices first.
- **Bit decomposition is not canonical.** 256 bits over a 255-bit field allow
  both k and k + p; for the k in steps (4)–(5) this gives a second valid r.
  Fix: enforce that the bits represent an integer below p.
- Prototype code: hash-to-curve is try-and-increment, nothing is constant
  time, and there has been no security review.

## Next steps

- Batch proving and batch verification as in §6.2, measured against Tables 4
  and 5 (one proof per dealer for n − 1 recipients).
- The two fixes above.
- Divisor-based exponentiation (§5.5 mentions about 546 constraints).

## Notes

Most tests, the Weierstrass gadget, the single-MSM verifier and `round1` were
written with the help of Claude (Anthropic); the commit history records which.
