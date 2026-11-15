# K0 report: exact oracles — GO (2026-11-13)

Gate (SKELETON §2-K0, ratified): two-way cross-checks ≤ 1e-10 on
full-rank AND rank-deficient anchors; CP/TP on every constructed Choi;
JRSWW sanity I ≥ −ln F_sq(ρ, R̄(ρ_{F1F2})) − 1e-12 on 500 seeded random
states. Finishing gate: `uv sync --locked && uv run pyrefly check &&
uv run pytest -q` → 0 type errors, 25 passed (1 informational cvxpy
warning, see Limitations).

## What exists

Package `petz_witness` (src layout, uv + pyrefly strict + pytest):
- `linops` — tripartition S⊗F1⊗F2, partial trace, Hermitian functions
  with ONE kernel threshold (CONF_KERNEL_EPS = 1e-13), seeded
  Ginibre states of chosen rank, depolarize.
- `entropies` — S_vN, H2, p2, p3, CMI, I′₂ (all nats; I′₂ never clamped).
- `petz` — R^P two independent ways (matrix powers vs Choi), twirled R̄
  two independent ways (closed-form kernel ω/sinh(ω) vs Gauss–Legendre
  quadrature of ∫β0(t)cos(ωt)dt), CPTP completion outside supp(ρ_F1),
  CP/TP assertion on the Choi object.
- `fidelity` — F_sq (locked convention), −ln F_sq.
- `anchors` — ρ⁻, ρ⁺ (Carbon v5, verified), GHZ (q/cl), W, product,
  Bell-with-trivial-F1 (fidelity-convention pin).
- `sdp` — CPTP-optimal fidelity of recovery (Watrous fidelity SDP +
  Choi TP constraint), default solver CLARABEL (0.9 s vs SCS 11.4 s on
  the anchor pair, and tighter: GHZ 0.500021 vs 0.500560; exact 0.5).

## Gates that went RED during development (M4 evidence, real not synthetic)

1. `partial_trace` off-by-one in the traced-axis pairing → 19/25 tests
   red (anchors, cross-checks, sanity). Fixed; the asymmetric-dimension
   test (2,4,2) now holds the orientation.
2. **Twirl kernel factor**: Carbon v5 states the kernel in the full-log
   convention (ω′/2)/sinh(ω′/2); transcribing against our half-log phase
   ω produced (ω/2)/sinh(ω/2). The independent quadrature path went RED
   (5 tests); the identity was re-pinned against scipy.quad by hand:
   true kernel ω/sinh(ω) at ω = ½ln(λ_aμ_y/(λ_bμ_x)). This is the
   designed catch of the two-path rule (T4/O2) — and a warning for K1:
   every Carbon formula gets convention-pinned before use.
3. Choi input-factor basis mismatch (eigen-labeled Choi vs computational
   application) found by derivation review before first run; the
   power-vs-Choi gate would have caught it.

Mutation calibration of the sanity gate (run, not claimed): root-fidelity
convention slip → RED 100/100; buggy half-log kernel → survivor for the
sanity gate (worst margin +0.079) but distinguished by the quadrature
cross-check gate, which owns kernel correctness (recorded in the test
docstring).

## Rule-evidence chain (substantial change, per skill Reporting)

- T4/O2 (independent oracles): R^P power-vs-Choi; R̄ closed-vs-quadrature;
  both ≤1e-10 / ≤1e-8 — caught bug #2. Unresolved: none.
- O22 (formula sources): each module header cites JRSWW/BSW/FR/Watrous/
  Berta–Tomamichel with the convention; Bell anchor tests the F_sq end.
- O1/R-CONFIG: CONF_KERNEL_EPS, CONF_QUAD_T_MAX/NODES, CONF_SANITY_*,
  CONF_SDP_MAX_TOTAL_DIM — single declarations with `#:` origins.
- O27: assert_cptp = runtime detection on the Choi object; stated residue:
  map QUALITY is owed to fidelity oracles/anchors, not to this check.
- M4: see RED list above; survivor reported with its distinguishing gate.
- M7: 500-state sanity does not prove the theorem (proved in source);
  gate role = bug detection; domain stated in the test.
- R-PINNED-BOUND: anchor expectations from hand derivations (analysis
  note) and closed forms; kernel values pinned against scipy.quad.
- M1: finishing gate green (above).

## Limitations (material)

- SDP solver reports `optimal_inaccurate` on some random states
  (informational warning left visible on purpose); accuracy is bounded
  by CONF_SDP_ATOL = 5e-4 assertions, adequate for a reference column,
  NOT for 1e-10-grade cross-checks.
- The 500-state sanity ensemble is 3-qubit (2,2,2); larger tripartitions
  are exercised only by the (2,4,2) cross-check case so far — K2 will
  widen this.
- `petz_apply_power` and the Choi path share `func_hermitian` (the
  kernel threshold), so a threshold-level bug would be common-mode;
  independence holds for everything above that one helper (declared
  blind spot).

## K0 verdict: GO — downstream gates may run (K1 next).
