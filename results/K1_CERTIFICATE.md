# K1 report: certificate theory — GO (2026-11-13)

Gate (SKELETON §2-K1, ratified): proofs closed; zero violations on 1e4
random spectra; δ_max^(3)-bound = 0 (≤1e-10) on flat spectra of every
rank 1 ≤ r ≤ d; continuous decay under depolarizing p → 0. Finishing
gate: pyrefly 0 errors, 39 passed.

## What exists

- theory/K1_CERTIFICATE.md — proofs, re-derived in project conventions
  (K0 lesson applied): Lemma 1 (identity + one-sided bound, two lines);
  Lemma 2 (k=2 ceiling exact and COMPLETE via the ≤2-level KKT argument;
  HT closed form pinned as an independent implementation, not assumed);
  Lemma 3 (k=3 sound closed-form bound via size-biased θ-split +
  Chebyshev — full derivation, including the C(θ)w sharpening term);
  Lemma 4 (≤3 levels at fixed p2, p3).
- src/petz_witness/certificate.py — delta_gap, delta_max2 (exact),
  delta_max2_ht_closed_form (pin), delta_max3_bound (sound witness
  ceiling), delta_max3_exact (best-effort slack reference, NaN-honest),
  witness_u2 / witness_u3 = I'2 + ceilings on SF1, F1F2.

## Numbers (own measurements)

- k=2 insufficiency (Lemma 2(iv), now a tested result):
  δ_max^(2)(1/2, 4) = 0.2796765 nats = 0.4035 bits — irreducible purity
  slack on a GHZ-grade marginal at ZERO noise; δ_max^(2) = 0 iff
  p2 ∈ {1, 1/d}.
- Continuity (depolarized classical GHZ, SF1 marginal, nats):
  p = 1e-2 → bound 4.09e-2 (actual δ 2.15e-2, ×1.9);
  p = 1e-4 → 7.05e-4 (×1.6); p = 1e-6 → 9.84e-6 (×1.5).
  Tighter than Carbon v5's quoted 3.9e-5 bits at p = 1e-6; near the
  flat manifold the bound is within 2× of the true gap — the regime
  that decides K2 recall.
- Slack vs best-effort exact δ_max^(3) at generic (p2, p3): ×2.6–3.2
  (probe set); the bound is loose far from flat, tight near flat.
  Measured, reported, NOT assumed away (O27 residue of Lemma 3).

## Gates RED during development (M4, real)

- The slack test itself caught a scipy fsolve pitfall: MINPACK returns
  ier=1 ("converged") AT A NON-ROOT when the Jacobian is singular (zero
  step passes the xtol criterion; reproduced minimally and recorded in
  the module). Fix: every root candidate is residual-verified against
  the moment equations (_ROOT_RESIDUAL_ATOL = 1e-10); the no-candidate
  outcome is NaN, not 0 (0 would falsely read "flat").

Mutations run (not claimed): σ² := 0 inside the bound → soundness RED
1699/2000 random spectra (the 1e4-gate owns this class); m=1-only
δ_max^(2) → 0/2000 violations (HT family IS the maximizer; the full
enumeration is completeness insurance per Lemma 2(i), kept because the
insurance is free); dropping C(θ)w → bound only loosens (one-sided gate
asymmetry, recorded).

## Rule-evidence chain

- T1: fsolve ier semantics established by probe, not memory — found the
  false-convergence behavior.
- T2: expected values (HT family max, two-level enumeration) computed by
  a dev probe before implementation; kept as tests.
- T4/O2: δ_max^(2) has two independent implementations (enumeration vs
  HT closed form), pinned equal; the bound is verified against exact δ
  of 1e4 independent spectra.
- O22: every lemma states its derivation inline in theory/; no inherited
  formula without a convention pin; Kong–Valiant NOT cited (flagged
  unverified — the k=3 minimality claim rests on our own Lemma 4).
- O27: Lemma 3's residue (bound-vs-exact slack) is measured and
  reported; delta_max3_exact explicitly labeled NOT-sound-as-witness.
- M4: see RED list; survivors reported with owning gates.
- M7: "zero violations on 1e4" domain stated (Dirichlet α∈{0.2,1,5},
  d ≤ 16, all ranks); finite corpus does not prove the inequality — the
  proof is Lemma 3; the corpus hunts implementation bugs.
- M1: finishing gate green.

## Limitations (material)

- delta_max3_bound slack ×2.6–3.2 away from the flat manifold — K2 must
  report recall for BOTH U^(2) and U^(3) so the slack's cost is visible.
- delta_max3_exact is best-effort (can return NaN; can underestimate);
  never feeds the witness.
- p3 estimator variance O(4^|M|) stated, not yet designed — K3 scope.
- Spectra corpus is classical (Dirichlet); quantum marginal CONSTRAINTS
  (which (p2,p3) pairs are jointly realizable by actual tripartite
  states) are not exploited — a possible tightening left to K2/paper.

## K1 verdict: GO — K2 (validity map) may run.
