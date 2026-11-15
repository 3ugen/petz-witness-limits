# K3 report: shot budgets — NO-GO executed (2026-11-13; v4 amendment 2026-11-14)

Verdict first: **the recorded NO-GO branch fires** ("if shadow/MLE
dominates — record it; that result feeds GO-B; valuable: 'below this m,
just reconstruct'"). Same-data shadow/MLE reconstruction dominates the
moment witness on shot count at every regime we could reach. The
witness's surviving delta is the soundness direction, not shot economy
(see K4_NOVELTY.md crossover table).

## Data

- results/k3_bench_m4.csv, k3_bench.csv (m=6) — v1 (blocked-functional
  estimator) full runs: recall 0% for BOTH methods everywhere; kept as
  the record of a diagnosed estimator flaw (below).
- results/k3_bench_m4_v3.csv — the record run: m=4 (dims 2,4,2),
  6 states (4 positives: ghz-cl/cluster at p ∈ {0, 1e-3}; 2 negatives:
  GHZ-q, W), budgets N ∈ {3e4, 1e5, 3e5}, 6 reps, v3 witness estimator
  + jackknife shadow/MLE baseline.

## Headline table (m = 4, record run)

| N | ours recall | baseline recall |
|---|---|---|
| 3e4 | 0% | 4.2% |
| 1e5 | 0% | 33.3% |
| 3e5 | 0% | 54.2% (probe outside grid: certifies ghz-cl-4 at UCB 0.04) |

UCB coverage of the exact value: ours 100%, baseline 99.1%. False
certifications of negatives: 0/36 both (on this family; the plug-in
has no guarantee — see the FoR-SDP probe below).

N*(80%): ours — not reached (extrapolated ~5e6 from the σ̂²-sensitivity
analysis); baseline — not reached in grid (extrapolated ~1e6 at m=4
from the recall curve). GO ratio (base/ours) < 1: **the gate's GO
clause is not met in any direction; the NO-GO branch is the result.**

## Why (measured mechanism chain)

1. v1 estimator flaw (ours, fixed): block-averaging the NONLINEAR
   witness functional carries a Jensen bias that does not shrink with N
   (~1.3 nats at exact 0). Replaced by full-dataset incomplete
   U-statistics + moment-box propagation (v2), then covariance-aware
   σ̂² propagation (v3).
2. The irreducible cost (ours): the ceiling δ_max^(3) enters through
   σ² = p₃ − p₂² with sensitivity ∂bound/∂σ² ≈ 1/(2θ²p₂²) ~ 8 at
   p₂ = ½ — certification at τ-headroom ~0.1 nats needs
   SE(σ̂²) ≲ 2e-3, i.e. N ~ 5e6 at m = 4 even with the v3 estimator
   (measured SE(σ̂²) = 0.024 at N = 3e5).
3. The baseline at m = 4 reconstructs a d = 16 state cheaply (the
   distinct-outcome space is only 81·16); its plug-in CMI with
   delete-block jackknife certifies from N ≈ 3e5.
4. At m = 6 BOTH fail at N ≤ 3e4 (v1 run; scale analysis: baseline
   needs ~4^m, ours needs 2^m moment variance × the ×8 sensitivity —
   the witness's break-even m lies where both budgets are already
   impractical for the systems in scope).
5. Soundness asymmetry (the surviving delta): plug-in FoR-SDP on
   shadow reconstructions UNDERSTATES infidelity on negatives
   (−ln 0.43 vs true 0.69 at N = 3e4) — errs toward false
   certification; our witness can only refuse.

## Scope notes / unrun arms (honest)

- Common-randomized-measurement estimators (2304.12292) and
  derandomized allocation: not run. In the single-tripartition setting
  every random product setting feeds all four marginals at once
  (downward closure), so allocation choices move constants, not the
  σ²-sensitivity wall; recorded as a limitation, not evidence.
- Many-pairs amortization axis: not run (the OR-clause of the gate was
  not needed — the main axis resolved the verdict; the axis remains
  future work for the paper's discussion).
- Negatives family is small (2 states × reps); "no false certs" for
  the baseline is an observation on this family, not a guarantee —
  exactly the point the soundness delta makes.

## v4 amendment (2026-11-14, post-Carbon-v6 — exact estimators)

Carbon v6 identified the v3 SE(σ̂²) = 0.024 as a Monte-Carlo floor from
the triple-subsample cap (verified: kernel var 1724 / 5e6 → 0.019) and
supplied exact complete U-statistics (Newton–Girard on accumulated
marginal shadows; verified against brute force to 1e-16). The v4 record
run (results/k3_bench_m4_v4.csv, m6_v4.csv) replaces the old wall:

| m=4, N | ours recall | baseline recall |
|---|---|---|
| 1e4 | 0% | 0% |
| 3e4 | 0% | 12.5% |
| 1e5 | 0% | 37.5% |
| 3e5 | 0% | 50% |
| 1e6 | **29.2%** | **100%** |

The witness DOES certify with exact estimators (first non-zero recall;
the v1/v3 zeros were estimator artifacts); N*(80%): baseline = 1e6,
ours > 1e6 (extrapolated ~2–3e6). Dominance factor ≈ 2–3×, not the
retracted ×16/"5e6" of the v3 text. Coverage: ours 100%, baseline
98.3%; false certs 0/60 both. The "×8 sensitivity" is replaced by the
exact Lemma-3 asymptotics B(σ²) ≈ 1.4·(σ²/p₂²)ln((d−1)p₂³/σ²)
(log-divergent slope, measured 26–62 over σ² ∈ [1e-5, 1e-2]).

m=6 (incl. the genuinely sub-tomographic N=300 < 3^6): recall 0% for
BOTH methods at every N ≤ 3e5 — together with the analytic
N*(m, τ=0.1) > 3^m for all m ≤ 16 (v6 §3c), certification in the
N_settings < 3^m regime is closed off for both pipelines. Regime
relabel accepted: the m=4 budgets are full-basis Pauli tomography
(N ≫ 3^4 = 81), not sub-tomographic; the paper says so.

Verdict unchanged: NO-GO branch stands (reconstruction dominates);
the witness's value is the soundness direction, now with honest
magnitudes. Run 2 (near-threshold negative false-cert rates) was lost
to an OOM kill before any output was flushed — pending rerun; until it
lands, the soundness-delta prose stays at "systematic downward plug-in
bias" strength.

## Rule-evidence chain

- T4: every estimate judged against exact oracles; estimator
  unbiasedness gates in tests/test_sampling.py.
- M4 (real REDs): v1 Jensen bias caught by the benchmark itself;
  code-index bug (2·basis+bit) caught by the probability-table pin;
  the Path+str crash cost one 20-min run (rerun, fixed).
- O7: subsample caps named (CONF_PAIR/TRIPLE_SUBSAMPLE_FULL).
- M7: all "0 false certifications" claims carry their finite domain.
- O10: estimator versions v1→v2→v3 registered with provenance in
  sampling.py (post-diagnosis, pre-record-run each time).

## K3 verdict: NO-GO branch executed → feeds GO-B.
Terminal letter per SKELETON §5 with K2 = boundary-GO, K3 = NO-GO,
K4 = GO-B: the paper is the GO-B limits/methods note under its
pre-named title.
