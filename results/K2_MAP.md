# K2 report: validity map — BOUNDARY outcome (2026-11-13)

Data: results/k2_map.csv — 9680 points: 10 seeds × 8 ensembles ×
(ratified noise grid + aux points 2e-3, 3e-3 for boundary localization).
Witnesses: U^(2), U^(3) (pre-registered), U^(3)± (registered after the
2-seed smoke diagnosed the mechanism, BEFORE this 10-seed run — see
theory/K1_CERTIFICATE.md Lemma 5 provenance note). Thresholds: ratified
[C5] set. Finishing gate: pyrefly 0 errors, 47 tests green.

## Headline numbers

**Soundness: 0 violations / 9680 points** (for all three witnesses) —
as the theorem requires; the map doubles as a K0 regression.

**Recall (pool = non-trivial MI ≥ 0.5·ln2 AND recoverable −lnF ≤ 0.02):**

| ensemble | p=0 | 1e-3 | 2e-3* | 3e-3* | 5e-3 | ≥1e-2 |
|---|---|---|---|---|---|---|
| (a) stabilizer-Markov, U^(3) | 64% | 43% | 30% | 6% | 0% | 0% |
| (a) stabilizer-Markov, U^(3)± | 69% | 59% | 48% | 36% | 8% | 0% |
| (h) QEC/repetition, U^(3)± | 100% | 100% | 100% | 100% | 100% | 0% |
| any ensemble, U^(2) | 0% | 0% | 0% | 0% | 0% | 0% |

(e) Gibbs and (f) HJPW: pools populated, recall 0 at ALL p including
p=0 — the structural false-refusal families, exactly as predicted
(non-flat recoverable states are refused by design). (b), (c), (g):
pools empty (no recoverable non-trivial members — expected). (d): pool
1–3 states, never certified.

**p_c (max p with recall ≥ 20%), per-seed spread 0.0% everywhere:**
(a)/U^(3): 2e-3; (a)/U^(3)±: 3e-3; (h)/U^(3)±: 5e-3.

## The mechanism (measured, not conjectured)

On noisy classical GHZ at p = 5e-3: CMI = 0.0033, −lnF = 0.0024, but
δ(SF1)+δ(F1F2) = 0.0746 while the DISCARDED pair δ(F1)+δ(SF1F2) =
0.0713. The one-sided witness pays the full noise-flatness defect of
two marginals and forfeits a near-perfect cancellation — that, not the
Lemma-3 slack (×1.8) and not I'2, is the dominant cost. U^(3)± recovers
the certified part of the cancellation (δ_min^(2) of the discarded
marginals) and shifts p_c from 2e-3 to 3e-3 (ensemble a).

## Gate verdict — honest reading

- **Ratified GO (U^(2)/U^(3), recall ≥ 20% at p ≥ 5e-3): FAILS.**
  U^(3) recall at 5e-3 is 0 on every ensemble.
- **Pre-registered NO-GO ("recall ≈ 0 at ALL realistic noise"): DOES
  NOT FIRE.** U^(3) certifies 43%/30% at p = 1e-3/2e-3; the certified
  region is non-empty and perfectly reproducible (spread 0%).
- **U^(3)± formally meets the GO clause via ensemble (h) at 5e-3** —
  but the (h) pool there is ONE distinct state (px=0, pz=1: the fully
  dephased GHZ₄, which is separable). Claiming GO on a single
  separable state would be gaming the gate (the exact failure mode
  Carbon v5 warned about in gate design). Not claimed.

**Verdict: BOUNDARY — neither pre-registered branch fires.** The map's
substantive content: a sound, reproducible certified region exists for
stabilizer-Markov states with boundary p_c ≈ 2–3e-3 per-qubit
depolarizing (3e-3 with the registered U^(3)± upgrade; the ratified
5e-3 target is missed by a factor ≈ 1.7). Per O10 this is an owner
decision point: the pre-registered gate pair did not anticipate a
non-empty-but-smaller-than-targeted certified region.

## Options for the owner (recorded, not chosen)

1. **Boundary-GO (recommended):** accept the measured boundary AS the
   K2 result; the paper is literally about validity boundaries, and
   "sound witness, certified region bounded by p_c ≈ 3e-3 per-qubit
   depolarizing on stabilizer-Markov states" is its honest headline.
   GO-A remains reachable via K3 (budget delta in the p < p_c regime).
   Honesty box carries: U^(3)± post-smoke registration; (h) thin pool.
2. **Strict NO-GO:** execute GO-B (limits note). Defensible but
   discards a non-empty, mechanism-understood certified region.
3. **Iterate the witness** (sound δ_min at fixed (p2,p3); quantum
   marginal constraints) under a NEW pre-registration round before any
   verdict. Costs time; the ×1.8 Lemma-3 slack and the remaining
   cancellation gap bound the possible gain to roughly p_c ≤ ~6e-3.

## Rule-evidence chain

- T4: every map quantity against K0 oracles; soundness = standing
  regression of the whole stack.
- O10: gate-gap surfaced and routed to the owner instead of silently
  reinterpreting the gate (this file §Gate verdict).
- M7: "0 violations / 9680" domain = the 8 named seeded ensembles ×
  named grid; no claim beyond it.
- M4 (ensemble layer): product-rotation perturbations were discovered
  USELESS (LU-invariance — CMI unchanged to 1e-14) and replaced by
  entangling 2-local rotations; the test
  test_entangling_perturbation_moves_cmi pins the fix.
- O27 residue priced: U^(2) dead on the whole map (k=2 insufficiency
  now an empirical map-wide fact, not only a lemma).

## Limitations

- Noise model is per-qubit depolarizing only (named); dephasing-
  dominant noise would flatten spectra differently — p_c under other
  channels is unmapped (paper scope note).
- Systems ≤ 5 qubits (d ≤ 32); p_c scaling with marginal dimension is
  visible in the (a) 4-vs-5-qubit split (U3± at 5e-3: GHZ-cl-4 passes
  0.0847 < τ, GHZ-cl-5 fails 0.108) but not systematically mapped.
- SDP reference column omitted from the sweep (CLARABEL per-point cost);
  available on demand for subsamples.
