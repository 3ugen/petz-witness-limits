# SEMANTICS — p3-petz-witness

Owner decisions 2026-11-13 (three forks resolved explicitly: v1 = simulation +
validity map, hardware deferred to a separate decision; money ceiling $0 until
an explicit owner decision; consumer decision = go/no-go on paper #2).
Preemption reading done before this file was settled:
notes/2026-11-13-p3-preemption-reading.md (niche verified free; LMW no-go
absorbed into the design).

Change 2026-11-13 (O10, after Carbon v5 skeleton review — REDESIGN
accepted, see notes/2026-11-13-carbon-v5-p3-skeleton-analysis.md): witness
redefined from (I′₂, four-marginal flatness certificate) to the one-sided
pair (I′₂, moment-ceiling δ_max^(k) on SF₁ and F₁F₂ only, k ∈ {2,3} with
p₃ = Tr ρ³ admitted — purities alone provably insufficient); oracle
hierarchy led by twirled Petz R̄; accuracy judged by additive
soundness/recall, not ratios. Consumer decision (§1) unchanged.

## 1. Consumer and decision

The consumer is the owner of the research track (RESEARCH.md) in subsequent
sessions. The decision made from v1 output: **go/no-go on paper #2** — a
methods note "validity boundaries of a cheap second-moment surrogate of Petz
recoverability" (GO-A/GO-B/NO-GO ladder with pre-registered gates, same
culture as p2). Secondary consumer — the owner's decoder-track positioning
(PLAN.md): results feed a conclusions section, but do not drive the decision.

*Retracted by:* a hardware run becoming part of v1 (then the consumer also
needs a QPU launch spec with budgets), or a machine consumer (pipeline reading
witness outputs programmatically) — then typed configuration export is
required, a textual note is not enough.

## 2. Admissible and refused

Can: exact simulation of small systems (3–10 qubits, S/F₁/F₂ partitions);
exact oracles — von Neumann CMI I(S:F₂|F₁), twirled-Petz R̄ (primary),
unrotated Petz R^P (secondary), and CPTP-SDP fidelity of recovery, all in
locked conventions (F_sq, nats); the surrogate I′₂ =
S₂(SF₁)+S₂(F₁F₂)−S₂(F₁)−S₂(SF₁F₂) with moment ceilings δ_max^(k)(SF₁),
δ_max^(k)(F₁F₂) for k ∈ {2,3} (p₂ and p₃ = Tr ρ³ from the same Pauli
bitstrings), and finite-shot one-sided confidence bounds in the
sub-tomographic setting regime; validity maps (soundness/recall vs
oracles) over state ensembles and noise; seed reproducibility.

Refusals that are NOT errors (contract features):
- the verdict "I′₂ tracks Petz fidelity only on the flat-spectrum manifold
  and the certifiable region is empty under realistic noise" is a valid
  scientific outcome that completes v1 successfully (NO-GO branch → limits
  paper, see §5);
- negative I′₂ values are displayed and flagged, never clamped (LMW: no
  universal sign);
- refusing to output a verdict when the flatness certificate fails is a
  feature: the witness is the PAIR (surrogate, certificate), and "certificate
  failed → no claim" is the designed behavior, not an error path.

Must refuse: claiming I′₂ is a Rényi generalization of CMI (BSW eq. 1.3
disqualification); any universal sign/bound claim for I′₂ (LMW no-go);
transferring Fawzi–Renner to second moments outside the certified regime;
claims about "recoverability in general" rather than Petz/rotated-Petz
(Sutter–Renner examples); claims at the HPTP/virtual level (VQMC Example 4:
not a function of moments); any hardware-behavior claim (noise models are
simulated, no QPU in v1); any spend of money.

## 3. Consequences of being wrong

A wrong GO → owner writes a paper whose central claim dies in external
review (reputational cost, wasted weeks). A wrong NO-GO → a free niche
(confirmed free on 2026-11-13) is abandoned while the Petz-plateau line heats
up — opportunity cost, partly irreversible (preemption risk by others).
Wrong numbers in the validity map → flawed hardware design in v2 and wrong
budgets. Scale: one owner. Detectability: high — exact simulation provides
independent oracles (CMI, Petz fidelity computed two ways: direct channel
application and fidelity formula); stabilizer states give closed-form limits
where surrogate = oracle exactly. Reversibility: full (local computation
only; $0 ceiling means no spent budget to lose).

## 4. Declared answers (pre-code; replaced by observation from sources)

- Owners: ONE (Eugen; agent sessions are the owner's tools).
- Exposure: outward — only result files (CSV/PNG/MD) inside the project's
  own results/ directory; no network, no secrets, no external input data, no
  timers. Randomness — present, injected via seed. (The VQMC benchmark
  states, if reused, enter as a vendored static file with provenance, not a
  network fetch.)
- Expected constructs: numeric formulas with sources (O22: BSW 1403.6102
  eq. 1.3/8.11; Fawzi–Renner 1410.0664; JRSWW 1509.07127; LMW 1212.0248;
  flatness bound |I′₂−I| ≤ Σ(H₀−H_∞)); semantic constants (O1/R-CONFIG:
  grids, δ-thresholds, budgets); bounded loops over ensembles/grids (O7);
  IO/randomness at the edges (O6); quantified claims with pre-registered
  gates (M7); Python package uv + pyrefly strict + pytest, kernels
  cross-checked to 1e-10 against an independent implementation (M1, p2
  precedent).

## 5. Out of scope, and what stops us

Out of scope for v1: hardware runs and hardware noise calibration (v2, only
after a GO gate and an explicit owner budget decision); realizing Petz
recovery as a circuit (that is the 2504.20399/2606.12020 line — we witness,
not implement); HPTP/virtual recovery certification (provably not
moment-accessible); von Neumann CMI estimation from shadows; Gaussian
continuous-variable branch (cited, not built); the old-project backlog (#5
replication, submit) — that stays in p2.

Direction stop (pre-registered before any results; reworded 2026-11-13
with the K2 gate after Carbon v5 — the multiplicative criterion was
ill-posed): if across the tested ensembles (SKELETON §2-K2 list a–h) the
witness U^(k) has recall ≈ 0 at every realistic noise level — i.e. the
certificate never fires on anything beyond trivial corners (product /
extreme-rank states) — then "cheap witness" is rejected as a construction,
and v1 concludes as the NO-GO limits paper: "Low-moment surrogates of
twirled-Petz recoverability: validity boundaries and shot-budget limits
against shadow reconstruction". (Soundness cannot fail — it is a theorem;
a soundness violation in the data is a code bug, not a scientific
outcome.) That outcome still satisfies the §1 consumer decision (a no-go
with evidence, not a failed project).
