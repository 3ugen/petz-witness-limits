# Paper #2 skeleton (living spine + KPI gates), v2

v1 created 2026-11-13 BEFORE any code or results; v2 same date — Carbon v5
skeleton review (reviews/google-carbon-v5/response.md) verdict REDESIGN
accepted after independent verification of every counterexample and a
citation audit (notes/2026-11-13-carbon-v5-p3-skeleton-analysis.md: math
clean, bibliography partly broken — no Carbon reference enters unverified).
Everything below remains prospective. Reread after every task; corrections
in §6 with a date. Paper draft language: English.

**Conventions (locked, K0 enforces):** fidelity is squared/Uhlmann
F_sq(ρ,σ) = ‖√ρ√σ‖₁²; logs are natural (nats) in all code and gates; bits
only in prose with explicit conversion. (Fawzi–Renner's F is root fidelity;
under F_sq the JRSWW bound reads I ≥ −ln F_sq(ρ, R̄(ρ_{F₁F₂})).)

## 1. Positioning

**Goal:** construct and calibrate a one-sided, low-moment upper bound on
twirled-Petz recovery infidelity,
U^(k) = I′₂ + δ_max^(k)(SF₁) + δ_max^(k)(F₁F₂),
where I′₂ = S₂(SF₁)+S₂(F₁F₂)−S₂(F₁)−S₂(SF₁F₂), and δ_max^(k)(M) is the
exact maximum of H₁(M)−H₂(M) consistent with the first k state moments
(k=2: purities, Harremoës–Topsøe extremal spectrum; k=3: adds p₃ = Tr ρ³
via the size-biased spectral variance σ² = p₃ − p₂², which vanishes iff
the spectrum is flat at any rank). Theorem chain (one-sided, unconditional):
−ln F_sq(ρ_{SF₁F₂}, R̄(ρ_{F₁F₂})) ≤ I(S:F₂|F₁) ≤ U^(k),
from the exact identity I − I′₂ = δ(SF₁)+δ(F₁F₂)−δ(F₁)−δ(SF₁F₂) with
δ(M) ≥ 0. All quantities estimable from single-qubit Pauli bitstrings in
the sub-tomographic shot regime (N_settings ≪ 3^|SF₁F₂|); map soundness
and recall against exact oracles (twirled Petz R̄ — primary; unrotated R^P
— empirical secondary; CPTP-SDP fidelity of recovery — reference) across
noisy stabilizer, Gibbs, HJPW, and W–GHZ ensembles; benchmark the
finite-shot upper confidence bound head-to-head against equal-budget
classical-shadow / projected-MLE CMI estimation from the same Pauli
dataset, including the many-pairs amortization axis (cost per certified
(F₁,F₂) pair over families of tripartitions).

**Working title:** "When can low-order Pauli moments witness twirled-Petz
recoverability? Validity boundaries and sub-tomographic measurement
budgets for a certified CMI surrogate" (revisit after K4).

**Honest framing (pre-committed):** validity-boundaries paper. Known
before any experiment: purities alone (k=2) provably cannot certify
intermediate-rank flat states (irreducible Θ(log d) slack, ≈0.404 bits
per GHZ-type marginal at zero noise) — the third moment p₃ is required,
so the pitch is "low-order moments", never "second moments suffice". The
witness is sound by theorem; the empirical questions are recall (how much
of the recoverable world it certifies) and shot cost.

**Two distinct terminal GOs (do not mix):**
- GO-A (positive methods paper): non-trivial certified recall under
  realistic noise (K2) + a sub-tomographic budget regime where the
  witness beats same-data shadow/MLE reconstruction (K3) + closed
  certificate theory (K1).
- GO-B (limits note): "Low-moment surrogates of twirled-Petz
  recoverability: validity boundaries and shot-budget limits against
  shadow reconstruction".

**Delta vs nearest works (preemption note + Carbon v5 §5, refs audited):**
- VQMC (2312.02031, IEEE TIT 2025): known-state SDP for virtual (HPTP)
  recoverability; we certify physical twirled-Petz recoverability from
  cheap data. Their Example 4 bars moment-based quantities from the HPTP
  level — our claims stay at Petz/CPTP. Mandatory positioning paragraph.
- Sutter–Renner (1705.06749): structural converse, nothing measurable;
  contrast formula D̄₂ with ω⁻¹ (their eq. 29).
- BSW (1403.6102) + LMW (1212.0248): the ceiling; I′₂ is NOT a Rényi-CMI;
  our ρ⁻/ρ⁺ anchors (classical Markov chains with I′₂ of both signs at
  F_Petz = 1) make the pathology concrete.
- Shadow/plug-in line: HKP (Nat. Phys. 16, 1050) reconstruction at small
  m; common randomized measurements (arXiv:2304.12292 = PRX Quantum 5,
  010352) — the honest same-data competitor family for K3, not an
  excluded topic.
- Decohered-stabilizer line (must cite; exact claims re-read by us before
  preprint): 2301.05238 (PRX Quantum 4, 030317 — Rényi-2 quantities on
  decohered toric code), 2310.08639 (PRX 14, 031044), 2404.07251
  (PRL 134, 070403 — Markov length, CMI decay, Petz-type quasi-local
  recovery; verified relevant), 2301.05687 (NOT yet verified by us).
- Fidelity of recovery: definitions Seshadreesan–Wilde 1410.1441;
  SDP/multiplicativity Berta–Tomamichel 1502.07973 (NOT 1412.4067, which
  is remainder terms — Carbon's ref corrected).
- Petz-implementation experiments (NMR 2025; ions 2504.20399; photonic
  2606.12020, tomography-verified): they implement, we witness.
- Torvinen et al. 2605.06848: Petz-plateau Darwinism — motivation for P8.

**NOT claimed:** universal witness (impossible from any finite moment set
— LMW); completeness (non-flat exact QMCs like ρ⁺ are refused — false
refusals are the designed price of soundness); unrotated-Petz theorem off
the flat manifold (open conjecture); HPTP/virtual level; "tomography-free"
as a blanket (the claim is sub-tomographic setting count, stated with the
regime); hardware (v1 simulation-only); "second moments alone".

**Target venue:** PRA / Quantum (methods).

## 2. Claims and KPI gates

Numeric thresholds below are owner engineering decisions; where a number
originated in Carbon v5 it is marked [C5], adopted provisionally 2026-11-13
and subject to owner ratification before the K2/K3 runs start.

### K0: oracles and conventions

Lock F_sq + nats. Implement CPTP-completed R^P and R̄ (kernel completion
Tr_{F₁}[(I−Π_supp)·]σ_S with explicit threshold 1e-13); R̄ via the
closed-form spectral kernel (ω/2)/sinh(ω/2) AND adaptive quadrature over
β₀(t) as an independent path; R^P via Kraus/channel application AND the
matrix-power formula. Exact oracles: von Neumann CMI; F_sq to R^P, R̄;
CPTP-SDP fidelity of recovery (Berta–Tomamichel) at 3–6 qubits. Anchors:
flat-spectrum states (surrogate exact; rotations cancel, R̄ = R^P — proven
only ON the manifold), ρ⁻/ρ⁺ (verified numbers in the analysis note),
W/GHZ, product states, rank-deficient marginals.
- GO: all two-way cross-checks ≤ 1e-10 on full-rank AND rank-deficient
  anchors; theorem sanity I ≥ −ln F_sq(ρ, R̄(ρ_{F₁F₂})) − 1e-12 on 500
  random states [C5].
- NO-GO: no downstream gate runs; fix or stop.
- **STATUS 2026-11-13: GO** — results/K0_ORACLES.md; 25 tests, pyrefly
  clean. Correction en route: the twirl kernel is ω/sinh(ω) in the
  half-log phase convention (Carbon v5 wrote the same kernel in the
  full-log convention; the transcription slip was caught by the
  independent quadrature path and re-pinned against scipy.quad). K1
  duty: convention-pin every inherited formula before use. SDP
  reference on CLARABEL; its accuracy class is 5e-4 (reference column
  only, not 1e-10 cross-check grade).

### K1: certificate theory

Prove: (i) the exact identity and the one-sided bound
I ≤ I′₂ + δ(SF₁) + δ(F₁F₂) (two-line, from δ ≥ 0); (ii) δ_max^(2)(p₂,d)
closed form on the Harremoës–Topsøe extremal spectrum (IEEE TIT 47,
2944–2960 (2001); Berry–Sanders J. Phys. A 36, 12255, quant-ph/0305059);
(iii) δ_max^(3)(p₂,p₃,d): exact maximum supported on ≤3 eigenvalue levels
(KKT on the simplex), machine-precision scalar search + the closed-form
continuous upper bound (θ-split + Chebyshev on the size-biased
distribution). Document the k=2 limitation as a result (0.404-bit GHZ
slack), not a footnote. p₃ estimator: 3-snapshot U-statistics on the same
Pauli data; variance O(4^|M|) — stated, with the shot-cost axis wired
into K2/K3. (Kong–Valiant-type k-moment slack: cite only if we verify the
exact theorem; currently flagged as extrapolation.)
- GO: proofs closed; zero violations on 1e4 random spectra;
  δ_max^(3) = 0 (≤1e-10) on flat spectra of every rank 1 ≤ r ≤ d;
  δ_max^(3) → 0 continuously as depolarizing p → 0 [C5].
- NO-GO: bound stays numerical observation; project slides toward GO-B.
- **STATUS 2026-11-13: GO** — results/K1_CERTIFICATE.md, proofs in
  theory/K1_CERTIFICATE.md (all re-derived, no inherited conventions;
  Kong–Valiant not cited — our own Lemma 4 carries k=3 minimality).
  Key measured facts: δ_max^(2)(½,4) = 0.2797 nats (insufficiency is a
  result); near-flat the k=3 bound is within ×1.5–1.9 of the true gap,
  far from flat ×2.6–3.2 (K2 must report recall for U^(2) AND U^(3)).
  New tool-trap recorded: scipy fsolve ier=1 at a non-root under
  singular Jacobian — every root is residual-verified now.

### K2: validity map (core gate — carries the direction stop)

Witness U^(k), k ∈ {2,3}. Judged by ADDITIVE tolerance, not ratio:
soundness (−ln F_sq(ρ,R̄) ≤ U^(k) everywhere — theorem; violations = code
bugs, wired to K0) and recall (share of truly recoverable states
certified). Pre-registered [C5, provisional]: noise grid
p ∈ {1e-3, 5e-3, 1e-2, 2e-2, 5e-2}; recovery threshold τ = 0.10 nats
(F_sq ≥ 0.905); triviality filter I(S:F₁) ≥ 0.5 ln 2; truly-recoverable
cut −ln F_sq ≤ 0.02.
Ensembles (seeded, named now): (a) stabilizer + coherent ε-perturbation;
(b) bounded-rank random / small Haar; (c) central-spin dephasing (p2
reuse); (d) LMW-adversarial corner; (e) 1D finite-T Gibbs chains
(XX+ZZ+field; CMI decay regime per Kato–Brandão CMP 370, 117 —
subexponential exp(−c√ℓ) — and Kuwahara–Kato–Brandão PRL 124, 220601,
above-threshold temperatures); (f) HJPW non-flat exact QMCs incl. ρ⁻, ρ⁺
and entangled block-diagonal with tunable q_j; (g) noisy W–GHZ
interpolations; (h) approximate-QEC states (S maximally entangled with a
small code under amplitude damping, F₁ = kept, F₂ = lost) — carries the
decoder-diagnostics positioning.
- GO: at p ≥ 5e-3 on ≥1 non-product entangled ensemble, recall ≥ 20%
  [C5] with zero soundness violations, and the certified-region noise
  boundary p_c stable to ≤10% across 10 independent ensemble seeds [C5].
  Report R^P and SDP fidelities alongside without claiming theorems.
- NO-GO (pre-registered direction stop = SEMANTICS §5): recall ≈ 0 at all
  realistic noise (certificate never fires on anything non-trivial) →
  execute GO-B; the map is the result. No third option.
- **STATUS 2026-11-13: BOUNDARY-GO (owner decision, option 1 of the
  report)** — the measured boundary IS the K2 result: sound witness,
  certified region bounded by p_c = 3e-3 per-qubit depolarizing
  (U^(3)±) on stabilizer-Markov states; headline carries the boundary;
  honesty box carries U^(3)± post-smoke registration and the rejected
  thin-pool (h) pass (results/K2_MAP.md). Soundness 0/9680; certified region
  non-empty, reproducible (p_c spread 0% across 10 seeds), boundary
  p_c = 2e-3 (U^(3)) / 3e-3 (U^(3)±, registered post-smoke pre-map) vs
  the ratified 5e-3 target; formal (h)-pass at 5e-3 rests on one
  separable state and is NOT claimed. Mechanism measured: the one-sided
  witness forfeits the δ(F1)+δ(SF1F2) cancellation under noise; U^(3)±
  recovers its certified part. U^(2): recall 0 map-wide. Options
  recorded in the report; recommendation: boundary-GO (option 1).

### K3: measurement-budget benchmark (sub-tomographic regime)

Known in advance (Carbon v5 blind-spot, accepted): deterministic Pauli
coverage of one tripartition's SF₁F₂ equals informationally complete
tomography (3^m settings) — so the benchmark lives at
N_settings ≪ 3^|SF₁F₂|, N_total ∈ [1e3, 1e5], m = |SF₁F₂| ∈ {3..10}, and
the estimator is the finite-shot upper confidence bound Û^(k)_{1−δ}
(δ = 0.05, shot variance of p̂₂, p̂₃ included) — point estimates fail
one-sided soundness ~50% at the boundary. Baselines on the EXACT same
Pauli bitstrings: shadow/projected-MLE reconstruction → plug-in CMI and
−ln F_sq(ρ̂, R̄(ρ̂)); common-randomized-measurement estimators
(2304.12292); random vs derandomized/importance Pauli allocation. Second
axis (our addition, p2 heritage): many-pairs amortization — cost per
certified (F₁,F₂) pair across overlapping fragment families, where
planned reuse is the structural advantage.
- GO: Û^(k)_{1−δ} keeps empirical soundness ≥ 1−δ and reaches 80%
  certified recall at m ≥ 6 with ≥1.5× fewer total shots than same-data
  shadow/MLE [C5], OR a measured ≥1.5× amortization delta on the
  many-pairs axis at equal per-pair error.
- NO-GO: if shadow/MLE dominates at all m ≤ 10 — record it; that result
  feeds GO-B (and is valuable: "below this m, just reconstruct").
- **STATUS 2026-11-13: NO-GO EXECUTED** (results/K3_BENCH.md). Shadow/
  MLE dominates at m=4 (certifies from N≈3e5; witness needs ~5e6 —
  σ²-sensitivity wall measured); at m=6 both fail ≤3e4 with scaling
  favoring reconstruction until impractical m. Estimator lineage
  v1→v2→v3 registered (Jensen-bias flaw diagnosed and fixed; the flaw
  was ours, the wall is the construction's). Surviving delta =
  guaranteed soundness direction (FoR-SDP-on-shadows measured erring
  toward false certification on negatives).

### K4: novelty = useful delta + citation hygiene

Quantitative crossover tables vs (1) VQMC SDP, (2) fidelity-of-recovery
SDP on shadow reconstructions, (3) shadow/plug-in CMI. GO: a concrete
(m, N_total, p) regime where the witness certifies with guaranteed
soundness while alternatives either lose soundness or need ≥1.5× shots;
otherwise GO-B. Citation duty before preprint: re-verify every reference
inherited from Carbon v5 (its bibliography was unreliable — see analysis
note), read 2301.05687, re-run the preemption search (Petz-decoder line
moves fast).
- **STATUS 2026-11-13: GO-B** (results/K4_NOVELTY.md). Crossover table
  built on K2/K3 data + the FoR-SDP unsoundness probe; no (m, N, p)
  regime beats alternatives on shots — the soundness delta and Lemmas
  1–5 + the validity map carry a limits note. Citation duties done:
  2301.05687 read (must-cite + companion 2301.05689; the p_c^(2)=0.178
  vs 0.109 window is a mandatory caveat); preemption re-search run —
  niche free; NEW mandatory refs: 2609.34699 (log-Euclidean Rényi CMI,
  Wilde group — theory window narrowing), 2610.12929 (Markov length
  instability). Terminal letter: GO-B under the pre-named title.

### K5 (placeholder, v2 of the project only): hardware

Unchanged from skeleton v1: out of scope; no gate text until an owner
budget decision — writing it now would fake pre-registration.

**Forbidden formulations (ledger-grade):** "I′₂ is a Rényi-2 CMI";
"I′₂ ≥ 0"; "small I′₂ ⇒ small CMI / recoverable" (universal); unrotated-
Petz theorem off the flat manifold; "Rényi-CMI monotone in α" without
caveats; "recoverability" without the twirled-Petz/CPTP qualifier;
"first"; HPTP/virtual-level claims; "second moments suffice" (k=2 slack
is a theorem); "tomography-free" without the sub-tomographic regime
qualifier; "classical F₁ restores monotonicity of I′₂" (refuted by ρ⁻);
multiplicative accuracy claims ("within factor x of −ln F").

**Allowed vocabulary:** "one-sided low-moment witness"; "sound by
construction, calibrated recall"; "coincides with I(S:F₂|F₁) on
flat-spectrum states"; "sub-tomographic setting budget"; "motivated by,
but distinct from, the Rényi-CMI of Berta–Seshadreesan–Wilde".

## 3. Paper structure

1. **Intro:** recoverability background (Petz, FR, JRSWW, SSA); cost of
   certifying recoverability today (tomography-verified Petz experiments;
   SDP routes); the question: what can low-order Pauli moments certify,
   at what setting/shot budget. Conventions stated up front.
2. **Ceiling and certificate (K1):** LMW no-go; exact identity; k=2
   closed form and its irreducible slack; k=3 variance certificate;
   accessibility costs.
3. **Validity map (K2):** soundness/recall over ensembles incl. Gibbs,
   HJPW (false-refusal anatomy via ρ⁻/ρ⁺), W–GHZ, QEC states.
4. **Budgets (K3):** sub-tomographic benchmark vs same-data shadow/MLE;
   many-pairs amortization; finite-shot one-sided CIs.
5. **Discussion:** decoder diagnostics (via ensemble (h)), Petz-plateau
   outlook (P8), HPTP boundary (what moments can never see); hardware as
   project-v2 outlook.

## 4. Task binding

- #1 → K0 (oracles, conventions, anchors incl. ρ⁻/ρ⁺, SDP reference)
- #2 → K1 (identity, δ_max^(2), δ_max^(3), estimator variances)
- #3 → K2 (ensembles a–h, map, direction-stop verdict)
- #4 → K3 (sub-tomographic benchmark + amortization axis)
- #5 → K4 (crossover tables, citation re-verification, preemption re-run)
- #6 → paper draft per §3, only after the K2 verdict is recorded

## 5. Reread protocol

After each task: update §2 statuses → reread §1 → on NO-GO execute the
recorded action → log in §6. Preprint gate: GO-A = K0 ∧ K1 ∧ K2(GO) ∧ K3 ∧
K4; GO-B = K0 ∧ K1 ∧ K2(NO-GO executed) ∧ K4. Owner ratifies [C5]
thresholds before #3 starts.

## 6. Correction log

- 2026-11-15 (v3.0): PAPER v0.2 — both draft reviews absorbed (Carbon
  v7 direct-API high: MAJOR REVISION; Qwen4-Max v7b budget-32k: MAJOR
  REVISION; convergent highest-value edit = claims-strength
  reconciliation). All checkable attacks verified before adoption
  (notes/2026-11-15-v7-draft-review-analysis.md): Θ(log d) scope,
  U = I on the §5.1 family (1e-14), Eq.(6) coefficient-1 asymptotics,
  slack decomposition ×8/×3.1, p_c definition reconciled with the
  recall curve, three-copy p₃ fix, 4^m rule removed, pre-registration
  sentence qualified, in-sample-coverage limitation stated. Figures
  1-2 generated from recorded CSVs; B.4's SE(σ̂²)@1e6 measured
  (3.7e-4) before printing. Remaining pre-submission duties: reference
  re-check, optional held-out coverage run, publication repo decision.

- 2026-11-14 (v2.9): #6 DONE — PAPER.md draft v0.1 written per §3
  structure under the adopted title; all numbers sourced from
  K0-K4 reports + Runs 1-3; honesty box (Appendix C) carries the five
  flagged deviations and the caught-bug list; references restricted to
  the audited set with a re-check note. Project tasks #1-#6 complete;
  remaining: owner read-through → external draft review (Carbon v7
  prompt on request) → arXiv/journal decision.

- 2026-11-14 (v2.8): Run 2 DONE on macOS arm64 (branch run2-mac-scs,
  merged; CSV re-verified here). Honest outcome: the anticipated
  decision-level unsoundness of plug-ins did NOT materialize — 0 false
  certs for all three pipelines; separation is in refusal margins
  (min UCB−τ: ours +0.45, MLE-CMI +0.19, FoR-SDP +0.05) and in the
  theorem-vs-observation status. Soundness-delta prose downgraded per
  the pre-committed rule (results/K3_BENCH.md §Run 2). Platform
  finding: CLARABEL SolverError on rank-deficient reconstructions
  (macOS arm64) — SCS fallback added to sdp.py. All review-mandated
  runs complete; #6 draft unblocked.

- 2026-11-14 (v2.7): Carbon v6 results review ACCEPTED (MORE RUNS
  NEEDED; analysis notes/2026-11-14-carbon-v6-results-analysis.md —
  5 numeric probes, incl. two confirmed defects of ours: min_floor
  grid soundness hole and the v3 Monte-Carlo SE floor). v4 exact-NG
  estimators implemented; Run 1 (K3 v4, m=4/m=6) DONE: witness
  certifies 29.2% at N=1e6 (m=4), baseline N*(80%)=1e6, dominance
  ≈2-3×; m=6 both 0% incl. N=300 — sub-tomographic regime closed for
  both. Run 3 (K2 amendments) DONE: p_c(τ) ≈ linear; p_c is
  channel-dependent — dephasing-biased noise extends it to ~1e-2.
  Run 2 (near-threshold negatives) LOST to OOM — rerun pending owner
  go (script hardened: incremental CSV, no tail-pipe). Title adopted:
  "One-sided low-moment witnesses of twirled-Petz recoverability:
  noise validity boundaries and shot-budget limits against Pauli
  state reconstruction". GO-B letter reconfirmed by the reviewer.

- 2026-11-13 (v2.6): #4/K3 DONE — NO-GO branch executed; #5/K4 DONE —
  GO-B (reports K3_BENCH.md, K4_NOVELTY.md). §1 reread: the spine
  ("validity boundaries … and sub-tomographic budgets") survives with
  the budget half answering NEGATIVELY — which the GO-B title already
  anticipated. Remaining task: #6 paper draft per §3 under the GO-B
  formula (K0 ∧ K1 ∧ K2 boundary-GO ∧ K3 NO-GO ∧ K4).

- 2026-11-13 (v2.5): owner decided K2 = boundary-GO (option 1). GO-A
  path continues via K3 in the p < p_c regime; #4/K3 started.
- 2026-11-13 (v2.4): #3/K2 RUN — verdict BOUNDARY, owner decision
  pending (results/K2_MAP.md; options recorded there). Design additions
  along the way, both registered with provenance: entangling 2-local
  perturbations replace useless product rotations in ensemble (a)
  (LU-invariance caught by test); witness U^(3)± = U^(3) −
  δ_min^(2)(F1) − δ_min^(2)(SF1F2) (Lemma 5, post-smoke pre-map). §1
  reread: the spine ("validity boundaries") is EXACTLY what the map
  produced — the open question is which GO letter carries it.
- 2026-11-13 (v2.3): #2/K1 DONE — GO (results/K1_CERTIFICATE.md). §1
  reread: spine holds; one sharpening — K2 will report recall for both
  U^(2) and U^(3) so the measured ×2.6–3.2 off-flat slack of the k=3
  bound is priced, not hidden. Witness functions witness_u2/u3 exist;
  K2 is unblocked.
- 2026-11-13 (v2.2): #1/K0 DONE — GO (results/K0_ORACLES.md). §1 reread:
  spine unchanged; one factual correction recorded in the K0 status
  (twirl kernel convention) strengthening the two-path rule for K1.
- 2026-11-13 (v2.1): owner ratified all [C5] thresholds as-is (K2: noise
  grid, τ = 0.10 nats, recall ≥ 20%, −ln F_sq ≤ 0.02 cut, p_c ≤ 10% over
  10 seeds; K3: δ = 0.05, ≥ 1.5×, 80% recall at m ≥ 6; K0: 500-state
  sanity; K1: 1e4 spectra). Task #1/K0 started.
- 2026-11-13 (v2): Carbon v5 REDESIGN accepted in full after independent
  verification (all counterexamples reproduced numerically; see analysis
  note). Changes: witness redefined as one-sided U^(k) on two marginals
  (four-marginal Σ(H₀−H∞) bound killed — discontinuous, wrong sign
  structure, inaccessible); conventions locked (F_sq, nats; DOC-3
  fidelity-convention bug fixed); K2 gate additive (soundness/recall),
  multiplicative factor-2 gate killed (ratio → 4 counterexample verified);
  R̄ promoted to primary oracle with closed-form kernel; k=3 moment p₃
  admitted (k=2 insufficiency is now a stated result); K3 moved to
  sub-tomographic regime with same-data shadow baselines + our
  amortization axis; ensembles e–h added; GO-B title fixed (forbidden-list
  violation); citation audit: Carbon's math clean, bibliography unreliable
  — refs corrected (1502.07973, 2304.12292, PRX 14 031044, CMP 370 117,
  quant-ph/0305059), unverified ones flagged. Thresholds marked [C5]
  pending owner ratification.
- 2026-11-13 (v1): skeleton created after SEMANTICS and preemption
  reading; all gates prospective; no results exist yet.
