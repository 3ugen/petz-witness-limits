# K4 novelty ledger (2026-11-13; preemption re-search same evening)

Verdict summary at bottom. Every reference here either verified today
(agent audit with full-text/abstract reads) or carried from the audited
2026-11-13 morning reading; nothing inherited unverified (Carbon v5
bibliography lesson applied).

## Crossover table: the witness vs its three nearest alternatives

Soundness = can the method falsely certify a non-recoverable state?
Our witness: NO by theorem (U^(3)± ≥ I ≥ −ln F_sq(ρ, R̄) always; 0
violations / 9680 K2 points + 0 false certs in K3 runs). Alternatives:

| Method | Access | Soundness | Shots (m=4, measured) | Regime limits |
|---|---|---|---|---|
| **U^(3)± UCB (ours)** | p₂×4, p₃×2 from random Pauli snapshots | exact at population level; asymptotic (1−δ) UCB on non-degenerate states | certifies from N ≈ 1–3e6 (29.2% recall at 1e6, v4 exact estimators) | p < p_c(τ=0.10, depol) ≈ 3e-3, up to ~1e-2 under dephasing-biased noise; m small |
| Shadow/MLE plug-in CMI + jackknife UCB (same data) | full reconstruction d=2^m | none (bias direction uncontrolled); empirically no false certs on our set | **N*(80%) = 1e6** (50% at 3e5; UCB 0.04 on the easiest state at 3e5) | needs ~4^m snapshots; m=6 fails ≤ 3e5 |
| FoR-SDP on shadow ρ̂ (Seshadreesan–Wilde def., Berta–Tomamichel SDP) | reconstruction + SDP | **unsound direction measured**: on the GHZ-q negative −ln F_sq = 0.43 (F̂_sq = 0.65) vs exact 0.693 (F_sq = 0.50) at N=3e4; 0.52 vs 0.693 at 1e5 — systematic downward plug-in bias on entangled negatives; Run 2 (near-threshold, 2026-11-14): 0 false certs for all pipelines at tested budgets; separation shows in refusal margins — min UCB−τ: ours +0.45, MLE +0.19, FoR-SDP +0.05 | similar to shadow/MLE | SDP cost caps m |
| VQMC SDP (2312.02031) | full classical description of ρ | different object (HPTP/virtual level; their Example 4: not moment-accessible) | not data-driven at all | — |

The honest K4 conclusion matching K3's data: at small m, straightforward
reconstruction DOMINATES on shot count; the witness's unique value is
the GUARANTEED soundness direction (and the K1/K2 theory: exact
ceilings, validity map, p_c boundary) — not shot economy. That is a
GO-B-grade delta: real, measurable (the FoR-SDP unsoundness probe),
but not the "≥1.5× cheaper" GO-A clause.

## Preemption status after the evening re-search (agent, 9 query lines)

Niche "sound witness from (p₂, p₃) + noise-validity map + budgets"
still free. New MANDATORY citations found:

- **2609.34699** Rubboli–Arqand–Wilde, "Log-Euclidean Rényi Conditional
  Mutual Information" (Sep 2026) — closest theory neighbor: Rényi-QCMI
  with Petz-recovery bounds; purely information-theoretic, the
  log-Euclidean form is NOT moment-accessible, no measurements/sample
  complexity. Explicit differentiation paragraph required. (Window is
  narrowing on the theory side — Wilde's group is building the
  Rényi-CMI apparatus; the measurement part should be secured fast.)
- **2610.12929** "An Instability of the Quantum Markov Length" (Oct
  2026) — CMI decay can switch exponential→algebraic under weak local
  channels; directly relevant caveat for our noise-validity map.
- **2610.11235** Markov length divergence in Markovian-physics systems
  — context for the map (nice-to-cite).
- **2301.05687** Bao–Fan–Vishwanath–Altman (read 2026-11-13): must-cite
  — canonical decohered-stabilizer line + the p_c^(2) = 0.178 vs
  p_c^(n→1) = 0.109 precedent. Direction corrected per v6: in the
  (0.109, 0.178) window the RAW Rényi-2 combination I′₂ ≈ 0 falsely
  signals decodability (unsound), while U^(3)± ≥ I forces refusal —
  the window is a MOTIVATION for the δ-ceilings, not a silence caveat.
  Companion **2301.05689** (Fan et al.: coherent information
  diagnostics) — must-cite alongside.
- Confirmed non-threats (implementation line, not certification):
  2610.20515 (Petz via single QSVT), 2511.05941 (Petz for optical
  losses), 2510.26895 (when is Petz cost-free), 2609.36080 (witness
  for process non-Markovianity — methodological parallel, cite),
  2611.10330 (adaptive Rényi-moment budgets — a TOOL for our K3 line,
  ally not rival), 2606.29443 / 2610.12320 (flatness-adjacent but
  other targets).
- 2605.06848 (Torvinen, Petz-plateau Darwinism) — the natural CONSUMER
  of the witness; mandatory in the discussion.

## What remains uncontested (the paper's deltas, GO-B grade)

1. A closed-form, sound-by-construction twirled-Petz recoverability
   certificate from second and third Pauli moments, with exact k=2 and
   closed-form k=3 ceilings (Lemmas 1–5, incl. the δ_min subtraction;
   uniqueness NOT claimed — moment-problem SDP bounds exist).
2. The quantitative validity map: non-empty certified region with a
   measured, reproducible noise boundary p_c ≈ 3e-3 — and the measured
   mechanism (forfeited δ-cancellation) behind it.
3. The budget reality check: head-to-head same-data comparison showing
   where reconstruction dominates and why (σ²-sensitivity), plus the
   measured unsoundness direction of plug-in FoR-SDP — a warning the
   community does not currently have in print.
4. k=2 insufficiency as a theorem + map-wide empirical fact.

## K4 verdict: GO-B (per the recorded gate: no (m, N, p) regime found
where the witness beats alternatives on shots; the soundness delta and
items 1–4 carry a limits/methods note). Final gate letter for the
project: GO-B with the pre-named title "Low-moment surrogates of
twirled-Petz recoverability: validity boundaries and shot-budget limits
against shadow reconstruction".
