# One-sided low-moment witnesses of twirled-Petz recoverability: noise validity boundaries and shot-budget limits against Pauli state reconstruction

**Draft v0.2 (2026-11-15). Target: Quantum (methods and limits); PRA
after APS formatting. v0.1 → v0.2 absorbs two independent adversarial
reviews (both MAJOR REVISION); every numerical claim below was
re-verified against the recorded data after the revision.**

## Abstract

Certifying that a tripartite state ρ_{SF₁F₂} is recoverable from
ρ_{F₁F₂} by the universal twirled Petz map on F₁ requires an upper
bound on the conditional mutual information I(S:F₂|F₁). Raw Rényi-2
entropy combinations cannot provide one: they can be negative or
vanish on non-recoverable states. We construct an exact one-sided
witness U ≥ I from second and third Pauli moments of marginals, by
combining the Rényi-2 combination with closed-form ceilings and floors
on the entropy gaps H₁ − H₂ at fixed moments. Purities alone provably
fail: their ceiling vanishes only at p₂ ∈ {1, 1/d} and grows as
(1 − r^{−1/2}) ln d for bounded-rank flat marginals; the third moment
repairs this through the size-biased variance σ² = p₃ − p₂², which
vanishes exactly on flat spectra. Across 9,680 exact-moment
state–noise evaluations we map the certified region and its
channel-dependent noise boundary, and trace the boundary to a
measured cancellation of entropy gaps that any one-sided moment bound
forfeits. On identical random-Pauli data, exact complete-U-statistic
certification loses to projected-MLE reconstruction by a factor 2–3
in shots at m = 4, and neither method certifies with fewer snapshots
than Pauli bases at m = 6: a pre-registered negative result
delimiting when low-degree moment witnesses can replace small-system
state reconstruction.

## 1. Introduction

Approximate quantum Markov chains sit at the junction of several
active lines: Petz recovery maps as near-optimal decoders, mixed-state
phases diagnosed by information quantities on decohered stabilizer
states [Lee–Jian–Xu, PRX Quantum 4, 030317 (2023); Sang–Zou–Hsieh,
PRX 14, 031044 (2024); Sang–Hsieh, PRL 134, 070403 (2025);
Bao–Fan–Vishwanath–Altman, arXiv:2301.05687; Fan–Bao–Altman–
Vishwanath, arXiv:2301.05689], and the operational reading of quantum
Darwinism as recoverability of a system from environment fragments
[Torvinen et al., arXiv:2605.06848]. The quantitative anchor is the
Fawzi–Renner theorem in its universal-recovery form [Fawzi–Renner,
CMP 340, 575 (2015); Junge–Renner–Sutter–Wilde–Winter, Ann. Henri
Poincaré 19, 2955 (2018)]: for the twirled Petz map R̄ built from
ρ_{SF₁} and ρ_{F₁},

    −ln F²(ρ_{SF₁F₂}, R̄_{F₁→SF₁}(ρ_{F₁F₂})) ≤ I(S:F₂|F₁).            (1)

(F² is squared/Uhlmann fidelity; all information quantities are in
nats; Fawzi–Renner's F is root fidelity, which accounts for the
absent factor 2.)

Certifying recoverability through (1) requires an upper bound on the
CMI, whose von Neumann entropies are not functions of any finite set
of polynomial moments. Petz-recovery experiments to date verify
recovery by full tomography [arXiv:2606.12020; arXiv:2504.20399]. At
the cheap end, the Rényi-2 combination

    I′₂ = H₂(SF₁) + H₂(F₁F₂) − H₂(F₁) − H₂(SF₁F₂)                    (2)

is emphatically NOT a conditional mutual information: Berta,
Seshadreesan and Wilde explicitly disqualify such linear combinations
as Rényi generalizations of CMI [JMP 56, 022205 (2015)], and
Linden–Mosonyi–Winter prove there are no non-trivial universal linear
inequalities among Rényi marginal entropies beyond non-negativity
[Proc. R. Soc. A 469, 20120737 (2013)]. We exhibit the pathology on
the simplest possible states: two classical Markov chains (perfectly
Petz-recoverable, CMI = 0) with I′₂ = ln(8/9) < 0 and
I′₂ = ln(10/9) > 0 (§3.1, Appendix A). In decohered stabilizer codes
the two-threshold picture of the mixed-state literature — Rényi-2
diagnostics seeing their own transition (≈0.178 for the toric code
under bit flips) while optimal decodability is lost at ≈0.109 [Fan et
al., arXiv:2301.05689; Sang–Zou–Hsieh, PRX 14, 031044] — warns that
quantities built from second moments alone can keep signaling
"Markovian" after true decodability is gone. An honest moment-based
recoverability statement must therefore be one-sided and must carry
its own correction terms.

Recent proper Rényi-CMI proposals — including the log-Euclidean Rényi
CMI of Rubboli, Arqand and Wilde [arXiv:2609.34699] — preserve
recovery bounds but involve non-commuting operator logarithms across
marginals, inaccessible from finite polynomial moments. Virtual
quantum Markov chain SDPs [Chen et al., IEEE TIT 71, 5387 (2025)]
characterize HPTP ("virtual") recovery from a full classical state
description; a CPTP certificate like ours trivially implies HPTP
recoverability, but neither CMI nor marginal moments can detect the
strictly wider virtual-recoverable regime where I(S:F₂|F₁) is large
(their Example 4 illustrates that virtual recoverability depends on
operator structure beyond CMI). The stability of CMI decay itself
under local noise is a live question [arXiv:2610.12929;
arXiv:2610.11235]. Here we ask the complementary measurement-level
question: **when can second and third Pauli moments provide a sound
one-sided certificate of physical twirled-Petz recoverability, which
noise mechanism sets the certificate's recall boundary even when the
true CMI remains near zero, and how do its finite-shot budgets
compare, on identical data, against simply reconstructing the state?**

Our contributions:

1. **A sound certificate, closed-form up to one scalar minimization**
   (§3): the exact identity I − I′₂ = δ(SF₁) + δ(F₁F₂) − δ(F₁) −
   δ(SF₁F₂), δ(M) = H₁(M) − H₂(M) ≥ 0, becomes I ≤ U after replacing
   the positive pair by moment ceilings and the negative pair by
   moment floors — the k = 2 extremizations are exact and complete
   (recovering the Harremoës–Topsøe and Berry–Sanders single-marginal
   forms), the k = 3 ceiling is a continuous closed-form bound in the
   size-biased variance σ² = p₃ − p₂². On exact stabilizer states all
   marginals are flat, every gap vanishes, and **U = I exactly**: the
   witness is tight precisely on stabilizer-Markov states.
2. **Purities alone provably do not suffice** (§3.2): the k = 2
   ceiling vanishes only at p₂ ∈ {1, 1/d}; for flat marginals of
   fixed rank r it grows as (1 − r^{−1/2}) ln d (0.2797 nats already
   at r = 2, d = 4), while for proportional rank r = cd it tends to
   −ln c. The viable certificate requires the third moment.
3. **The validity map** (§4): 9,680 exact-moment evaluations (121
   state instances × 8 noise levels × 10 seeds, eight ensembles). The
   certified region is non-empty; with p_c defined as the largest
   grid noise rate with ensemble recall ≥ 20%, the boundary on the
   stabilizer-Markov ensemble is p_c = 3×10⁻³ per-qubit depolarizing
   (m = 4–5, τ = 0.10 nats), extending to 10⁻² under 10:1
   dephasing-biased noise. On a classical-GHZ member we measure the
   mechanism directly: the forfeited cancellation of entropy gaps
   inflates the certified target ×8 above the true CMI, and the
   ceiling's slack a further ×3.
4. **A pre-registered negative budget verdict** (§5): with exact
   complete U-statistics, witness certification at m = 4 first
   appears at N = 10⁶ snapshots (29.2% recall) while projected-MLE
   reconstruction reaches 50% at 3×10⁵ and 100% at 10⁶; the measured
   dominance is a factor 2–3 wherever either method works (witness
   N*(80%) ≈ 2–3×10⁶ by extrapolation, Appendix B.4). At m = 6
   neither method certifies anything up to 3×10⁵ snapshots, including
   N = 300 < 3⁶.
5. **What survives** (§5.1): at exact moments, soundness I ≤ U is a
   theorem, and on our near-threshold negative family U = I exactly;
   at finite shots all three tested pipelines (witness UCB, MLE-CMI
   jackknife, FoR-SDP jackknife) produced zero false certifications,
   with refusal margins explained quantitatively by their different
   population estimands plus finite-shot width.

All gates, tolerances and fallback branches were fixed before the
runs they judge, with recorded deviations (Appendix C); the data,
seeds and scripts regenerate every number in this paper.

## 2. Setup and conventions

Tripartition S ⊗ F₁ ⊗ F₂, dimensions (d_S, d₁, d₂); all numerics use
qubit registers with contiguous blocks; m = |SF₁F₂| denotes the total
qubit count. Fidelity is squared/Uhlmann, F²(ρ,σ) = ‖√ρ√σ‖₁²; logs
natural. H_α is Rényi entropy, H₁ von Neumann, H₂ = −ln p₂ with
p₂ = Tr ρ², p₃ = Tr ρ³.

Twirled Petz map: R̄(·) = ∫ β₀(t) R^{t/2}(·) dt, β₀(t) =
(π/2)(cosh πt + 1)^{−1}, with rotated maps as in JRSWW. In the
eigenbases ρ_{SF₁} = Σ λ_a |a⟩⟨a|, ρ_{F₁} = Σ μ_x |x⟩⟨x|:

    R̄(Y) = Σ_{a,b,x,y} (ω/sinh ω) √(λ_a λ_b / μ_x μ_y)
            |a⟩⟨a| (I_S ⊗ |x⟩⟨x| Y |y⟩⟨y|) |b⟩⟨b|,
    ω = ½ ln(λ_a μ_y / (λ_b μ_x)),                                   (R̄)

completed to CPTP on (supp ρ_{F₁})^⊥ by routing kernel weight to
ρ_{SF₁}. The kernel ω/sinh ω at half-log phase was verified against
direct quadrature of β₀ (a half-log/full-log convention slip is easy
to make and was caught in our own pipeline by this two-path check).
All oracles are implemented twice through independent code paths
agreeing to 10⁻¹⁰.

Noise channels (Kraus conventions used throughout): per-qubit
**replacement depolarizing** E_p(ρ) = (1−p)ρ + p I/2 per qubit
(equivalently Kraus {√(1−3p/4) I, √(p/4) X, √(p/4) Y, √(p/4) Z});
**10:1 dephasing-biased Pauli** with p_Z = 10 p_X = 10 p_Y and total
error probability p; **amplitude damping** with per-qubit γ = p. When
comparing p_c values to the stabilizer-code literature (which uses
the Pauli channel (1−p)ρ + (p/3)(XρX+YρY+ZρZ)), rates convert by a
factor 3/4 in bit-flip probability.

## 3. The certificate

### 3.1 Identity and one-sided bound

**Lemma 1.** For any ρ on S ⊗ F₁ ⊗ F₂, with δ(M) := H₁(M) − H₂(M) ≥ 0:

    I(S:F₂|F₁) − I′₂ = δ(SF₁) + δ(F₁F₂) − δ(F₁) − δ(SF₁F₂)           (3)

exactly; hence, with (1),

    −ln F²(ρ, R̄(ρ_{F₁F₂})) ≤ I ≤ I′₂ + δ(SF₁) + δ(F₁F₂).            (4)

Proof: both sides of (3) are the same ±combination of marginal
entropies at α = 1 and α = 2; δ ≥ 0 is monotonicity of α ↦ H_α.
Only SF₁ and F₁F₂ need upper certificates; the discarded
−δ(F₁) − δ(SF₁F₂) ≤ 0 strengthens (4) — or is partially reclaimed
(Lemma 4). The classical anchors ρ⁻, ρ⁺ (Appendix A) show I′₂ is
sign-indefinite on perfectly recoverable states.

### 3.2 Moment ceilings and floors

**Lemma 2 (k = 2 ceiling, exact and complete).** The maximizer of H₁
at fixed (p₂, d) has at most two distinct nonzero eigenvalues.
*Proof sketch* (Appendix B.1): on the support, Lagrange stationarity
forces ln λ to be affine in λ, and a strictly concave function meets
an affine one in ≤ 2 points; zero eigenvalues recurse to lower rank;
the only configuration with linearly dependent constraint gradients
is the uniform spectrum, enumerated separately. The maximum is the
one-high-rest-flat family a = (1+√((d−1)(dp₂−1)))/d [Harremoës–
Topsøe, IEEE TIT 47, 2944 (2001); Berry–Sanders, J. Phys. A 36,
12255 (2003)]; we write δ̄⁽²⁾(p₂, d).

**Proposition 1 (k = 2 insufficiency).** δ̄⁽²⁾(p₂, d) = 0 iff p₂ = 1
or p₂ = 1/d. For flat states of fixed rank r (p₂ = 1/r),
δ̄⁽²⁾(1/r, d) ~ (1 − r^{−1/2}) ln d as d → ∞ (already 0.2797 nats =
0.4035 bits at r = 2, d = 4), while at proportional rank r = cd,
δ̄⁽²⁾ → −ln c (e.g. ln 2 at r = d/2), and at r = d − 1 it vanishes as
O(1/d). Purity-only certification of a GHZ-type (bounded-rank)
marginal is impossible at any noise level; the third moment is
required.

**Lemma 3 (k = 3 ceiling, closed-form and sound).** With
σ² := p₃ − p₂² ≥ 0 (variance of the size-biased spectral variable;
σ² = 0 iff the spectrum is flat at some rank), for every θ ∈ (0,1):

    δ ≤ σ²/(2θ²p₂²) + max_{0≤w≤η} [ w ln((d−1)p₂/w) − C(θ) w ],      (5)

η = min(1, σ²/((1−θ)²p₂²)), C(θ) = (1−θ) + (1−θ)²/(2θ²) (proof:
Appendix B.2; the inner maximum is explicit, g(w*) = w* at
w* = (d−1)p₂ e^{−1−C}). The bound vanishes on every flat spectrum of
any rank and is continuous in (p₂, p₃). We write δ̄⁽³⁾(p₂, p₃, d) for
the θ-minimized value (a one-dimensional minimization of an explicit
function; any fixed θ already gives a valid bound). Its small-σ²
behavior is

    δ̄⁽³⁾ ~ (σ²/p₂²) · L,   L := ln((d−1)p₂³/σ²),  σ² → 0,           (6)

with limiting coefficient 1; across the window relevant in §5
(L ≈ 4.5–11) the measured ratio to the leading term is 1.36–1.47, and
we use 1.4 as a conservative envelope at (p₂ = ½, d = 8). The
log-divergent slope ∂δ̄⁽³⁾/∂σ² controls the shot budgets of §5.

**Lemma 4 (k = 2 floor, exact).** The same two-level enumeration
computes δ̲⁽²⁾(p₂, d) = min δ at fixed p₂, attained at m₀ = ⌊1/p₂⌋
copies of one level plus a single smaller one, and vanishing at every
p₂ = 1/r. Since the minimum at fixed (p₂, p₃) is no smaller,

    I ≤ U := I′₂ + δ̄⁽³⁾(SF₁) + δ̄⁽³⁾(F₁F₂) − δ̲⁽²⁾(F₁) − δ̲⁽²⁾(SF₁F₂) (7)

remains sound. U is the witness throughout; "certify" means U ≤ τ
with the pre-registered τ = 0.10 nats (F² ≥ 0.905 for R̄).

On exact stabilizer states all four marginals are flat, so every
ceiling and floor vanishes and U = I exactly — the certificate is
TIGHT on the stabilizer-Markov manifold, and all its slack is a
noise-robustness price. In the decohered-code two-threshold window
mentioned in §1 we expect the ceilings to force refusal where raw
I′₂ ≈ 0 would mislead; we have not evaluated toric-code marginals
here, and note that δ(M) grows with |M| under local noise, so U is a
small-subsystem tool, not a phase-boundary diagnostic.

## 4. Validity map

Eight seeded ensembles (Appendix A): (a) stabilizer-Markov bases with
entangling perturbations; (b) bounded-rank random; (c) central-spin
dephasing; (d) classical branching with pointer bias; (e) 1D Gibbs
chains; (f) non-flat exact Markov chains of Hayden–Jozsa–Petz–Winter
form [CMP 246, 359 (2004)]; (g) W–GHZ interpolations; (h) a
repetition-code family — 121 state instances × 8 noise levels × 10
seeds = 9,680 exact-moment evaluations on 4–5-qubit tripartitions.
Numerical regression: 0 violations of I ≤ U across all points (as
the theorem demands; this cross-checks the independent oracle
implementations, and is not an empirical finding about states).

Definitions: the recall pool is recoverable (−ln F² ≤ 0.02) and
non-trivially correlated (I(S:F₁) ≥ ½ ln 2) states; **p_c := the
largest grid noise rate with pool recall ≥ 20%** (grid {0, 10⁻³,
2×10⁻³, 3×10⁻³, 5×10⁻³, 10⁻², ...}).

Results on the stabilizer-Markov ensemble (m = 4–5, τ = 0.10,
replacement depolarizing; Fig. 1):

| p | 0 | 10⁻³ | 2×10⁻³ | 3×10⁻³ | 5×10⁻³ | 10⁻² |
|---|---|---|---|---|---|---|
| U⁽²⁾ (purities only) | 0 | 0 | 0 | 0 | 0 | 0 |
| U⁽³⁾ (no floors) | 64% | 43% | 30% | 6% | 0 | 0 |
| U (Eq. 7) | 69% | 59% | 48% | 36% | 8% | 0 |

- **p_c(τ=0.10, m=4–5, depolarizing) = 3×10⁻³**, with the identical
  grid crossing in all 10 seeds (i.e. reproducible within the grid
  spacing [3×10⁻³, 5×10⁻³)). Under 10:1 dephasing-biased noise the
  certified region extends to **p_c = 10⁻²** (recall 35% at 5×10⁻³,
  21% at 10⁻²); amplitude damping tracks depolarizing (p_c = 3×10⁻³).
  The τ-dependence is monotone and roughly affine over τ ∈ [0.05,
  0.25] (p_c: 1.1×10⁻³ → 10⁻²).
- The purity-only witness certifies nothing anywhere — Proposition 1
  is a map-wide empirical fact as well as a theorem.
- At p = 0, every unperturbed stabilizer-Markov base state is
  certified with U = I = 0 exactly (the §3.2 tightness); the 69%
  headline reflects the entangling-perturbed members.
- **Mechanism, measured on the classical-GHZ member at p = 5×10⁻³**
  (consistent with, though measured on one member of, the ensemble
  behind p_c): true CMI 0.0033, recovery infidelity 0.0024, yet
  δ(SF₁)+δ(F₁F₂) = 0.0746 against a discarded cancellation
  δ(F₁)+δ(SF₁F₂) = 0.0713. After the Lemma-4 floors reclaim part of
  it, an exact-gap bound would read 0.0273 (×8 above the true CMI),
  and the Lemma-3 ceiling adds ×3.1, giving U = 0.0847 — still
  certified on this member; the ensemble recall at 5×10⁻³ (8%) is
  lost on the larger and perturbed members. This cost is intrinsic
  to any one-sided bound built on (3): it bites while the true CMI
  is still negligible.
- Structural false refusals: Gibbs chains and non-flat HJPW exact
  Markov chains are recoverable but never certified (recall 0 even
  at p = 0). Soundness is preserved; completeness is the designed
  sacrifice — the witness detects flat-spectrum Markovianity, and we
  claim nothing broader.

## 5. Finite-shot budgets: witness vs reconstruction

Estimation uses random-Pauli classical shadows [Huang–Kueng–
Preskill, Nat. Phys. 16, 1050 (2020)]: each snapshot is one shot in
an i.i.d. uniform product Pauli basis. Both pipelines consume the
SAME bitstrings.

**Witness pipeline.** Exact complete U-statistics for p₂ and p₃ of
the four marginals via a Newton–Girard reduction on the accumulated
marginal shadow (Appendix B.3; cost O(min(N, 6^{|M|}) d_M²), no
subsampling), exact first-order influence functions for variances and
the p̂₂–p̂₃ covariance, the bias correction σ̂² = p̂₃ − p̂₂² +
V̂ar(p̂₂), and a Bonferroni moment box evaluated by worst-case search
over each marginal's confidence interval (including every flat point
1/r, where the floors vanish) — yielding a one-sided 95% upper
confidence bound Û. The UCB is asymptotically (1−α) calibrated on
non-degenerate states; on maximally mixed marginals the first-order
projection vanishes and the normal quantile is approximate — a caveat
that applies exactly on the stabilizer-like target states. Empirical
coverage of the exact U was 100% across all runs below; because the
estimator went through four recorded revisions on development data
(Appendix C), we additionally ran a held-out check on fresh seeds
never used in development (seven state families including the
degenerate maximally-mixed-F₁ case, two budgets, 350 runs):
coverage 350/350.

**Reconstruction pipeline.** Projected MLE (PSD, unit trace) of the
full ρ_{SF₁F₂} from the same shadows, plug-in CMI, delete-block
jackknife (20 blocks) 95% UCB. At m = 4 every budget below puts
≥ 120 shots into each of the 3⁴ = 81 Pauli bases: this is full-basis
Pauli state reconstruction, and we label it as such.

**Result (m = 4, dims 2×4×2; 4 recoverable positives and 2
negatives, 6 repetitions each — 24 positive and 12 negative runs per
budget; states and their exact values in Appendix A; Fig. 2):**

| N snapshots | witness recall | reconstruction recall |
|---|---|---|
| 10⁴ | 0 | 0 |
| 3×10⁴ | 0 | 12.5% |
| 10⁵ | 0 | 37.5% |
| 3×10⁵ | 0 | 50% |
| 10⁶ | 29.2% (7/24) | 100% (24/24) |

Defining N* as the smallest grid budget with ≥ 80% positive recall:
reconstruction crosses at N* = 10⁶ (recall 100% there, 50% at
3×10⁵); the witness first certifies anything at 10⁶ and its N* is an
extrapolation, ≈ 2–3×10⁶ from the variance model of Appendix B.4.
UCB coverage of the exact value: 100% (witness, 180/180), 98.3%
(reconstruction, 177/180); zero false certifications of negatives
for both (0/60 each over all budgets).

At m = 6, neither method certifies anything up to N = 3×10⁵,
including the genuinely sub-tomographic point N = 300 < 3⁶ = 729
total snapshots; the variance model of Appendix B.4 keeps the
witness's N* above 3^m through m ≤ 16 at τ = 0.10 (a model-based
statement with constants displayed there, not a measurement).

The driver is (6): with slope ∂δ̄⁽³⁾/∂σ² ≈ 1.4 L/p₂² ≈ 38 at
(p₂ = ½, d = 8, σ² ≈ 10⁻³), two ceiling marginals and the Bonferroni
quantile z ≈ 2.6 require SE(σ̂²) ≲ τ/(2·38·2.6) ≈ 5×10⁻⁴ — reached
only near N ~ 10⁶ at m = 4 — while reconstructing a d = 16 state
from ≥ 10⁵ full-basis shots is comparatively easy. **At m = 4, and
for every system size where projected MLE remains computationally
tractable, reconstruction dominates on single-copy Pauli shots by
the measured factor 2–3.** We looked for a regime where the witness
is cheaper and did not find one; this negative is a primary,
pre-registered result. The closest prior art on randomized-
measurement moment estimation [Vermersch et al., PRX Quantum 5,
010352 (2024)] optimizes estimators of the moments themselves; our
negative result is about the certification functional built on them
— the log-divergent ceiling sensitivity, not the moment estimators,
is the binding constraint.

### 5.1 Near-threshold negatives: what the guarantee structure buys

On the mixture family ρ(λ) = (1−λ)|GHZ₄⟩⟨GHZ₄| + λ ρ_GHZ-cl tuned to
exact twirled-Petz infidelity 0.13 and 0.18 (true negatives just
above τ = 0.10), all proper marginals are λ-independent and flat, so
the ceilings and floors vanish and **U(ρ(λ)) = I(ρ(λ)) exactly** —
the population margins are U − τ = I − τ = +0.133 / +0.209 nats,
versus −ln F² − τ = +0.030 / +0.080 for the fidelity-of-recovery
estimand [Seshadreesan–Wilde, PRA 92, 042321 (2015); SDP form after
Berta–Tomamichel, IEEE TIT 62, 1758 (2016)]. All pipelines bound
their statistic in nats; they differ in the estimand (CMI for the
witness and MLE pipeline, −ln F² for FoR-SDP), and (1) orders the
two.

Protocol: 100 repetitions per budget N ∈ {3×10³, 10⁴, 3×10⁴} for the
witness and MLE-CMI; the FoR-SDP arm ran 30 repetitions per target at
N = 10⁴ with an 8-block jackknife (SCS solver; cost-limited, stated
as such). Results: **zero false certifications for all three
pipelines** (witness 0/600, MLE-CMI 0/600, FoR-SDP 0/60). Minimum
UCBs over all runs: witness 0.548/0.636, MLE-CMI 0.294/0.312,
FoR-SDP 0.155/0.179 — each pipeline's margin is its population
margin plus its finite-shot width, and at these budgets (where the
witness certifies no positives) the witness's larger margin reflects
estimator width, not a population advantage: on this family the
population certificate margin U − I is exactly zero. A
far-from-threshold probe does show the known downward plug-in bias
of FoR-SDP on an entangled negative (estimated infidelity 0.43 vs
true 0.693 at N = 3×10⁴ on pure GHZ₄).

The honest summary: at exact moments the witness's soundness is a
theorem, and nothing in our finite-shot experiments contradicted its
asymptotic 95% calibration; the plug-in pipelines carried no
guarantee but were never caught false-certifying either. The
guarantee structure is a difference in what can be PROVEN, not — at
these budgets and on this family — a measured difference in decision
errors.

## 6. Discussion

**When to use what.** For m ≤ 6 and the budgets measured here,
reconstruction with a jackknifed plug-in reaches reliable
certification decisions at 3×10⁵–10⁶ random-Pauli snapshots (m = 4)
and is the pragmatic choice; the moment witness needs ~2–3× more and
offers in exchange an exact population-level inequality and exactness
on stabilizer-Markov states. Use the witness where (i) the one-sided
guarantee at exact moments matters more than shot economy — e.g.
certifying recoverability of tensor-network or otherwise classically
represented marginals, where moments are exact and the finite-shot
question is absent; (ii) the states are near the flat-spectrum
manifold at per-qubit error ≲ 3×10⁻³ (≲ 10⁻² under dephasing-biased
noise); or (iii) as an auditable screening layer whose refusals are
free of reconstruction bias. If the target state is known, direct
fidelity estimation [Flammia–Liu, PRL 106, 230501 (2011)] is the
cheap alternative to both. The by-product lemmas — exact k = 2
entropy-gap extremizations and the σ²-flatness characterization —
are independently useful wherever "is this marginal flat?" is the
real question.

**Open directions.** (i) Three-copy access: p₃ is a three-copy
observable; coherent cyclic-shift measurements would bypass the
single-copy shadow variance wall that sets our negative verdict, and
the known single-copy vs multi-copy separations for purity-type
functionals suggest this is where a moment witness could win.
(ii) Many-pairs amortization: across families of overlapping
(F₁, F₂) tripartitions the same snapshots feed every marginal; not
benchmarked here. (iii) Joint quantum constraints on the moment box
(the Bonferroni box treats marginals as independent). (iv) Mapping
p_c for witness classes beyond flat-spectrum, in light of the
instability of CMI decay under local channels [arXiv:2610.12929].

## Appendix A: anchors, ensembles, and benchmark states

Anchors: ρ⁻ = ½|00⟩⟨00|_{SF₁} ⊗ I/2 + ½ (I/2)_S ⊗ |10⟩⟨10|_{F₁F₂}
and ρ⁺ = ½|000⟩⟨000| + ½ (I/2)_S ⊗ |1⟩⟨1|_{F₁} ⊗ (I/2)_{F₂}:
classical Markov chains, CMI = 0, Petz fidelity 1, I′₂ = ln(8/9) and
ln(10/9).

Ensembles (per seed; dims (2,4,2) and (2,4,4)): (a) classical GHZ
and 1D cluster bases (4–5 qubits) alone and under entangling
two-local perturbations exp(−iεG) on adjacent pairs, ε ∈ {0.05, 0.1,
0.2}, 3 draws each (40 states); (b) Ginibre states of rank {1, 2,
d/2, d} (24); (c) central-spin dephasing, seeded couplings, 3 times
(9); (d) classical branching with pointer bias (7); (e) Gibbs chains
of H = Σ(J_x XX + J_z ZZ + h Z), β ∈ {0.5, 1, 2, 4} (16); (f) HJPW
block states Σ_j q_j ρ_S^j ⊗ |j⟩⟨j| ⊗ σ^j (6); (g) W–GHZ
interpolations (10); (h) reference ⊗ 3-qubit repetition code under
bit-flip/dephasing mixes (9). Total 121 instances; noise grid of 8
points; 10 seeds.

§5 benchmark states (m = 4): positives ghz-cl-4 and cluster-4 at
p ∈ {0, 10⁻³} (exact −ln F² ≤ 0.0005, exact U ≤ 0.037); negatives
|GHZ₄⟩ (I = ln 2, −ln F² = 0.693) and |W₄⟩ (−ln F² = 0.284). §5.1
family as in the text with λ tuned by bisection to the stated exact
infidelities (λ = 0.346, 0.258).

## Appendix B: proofs, estimators, and the variance model

**B.1 (Lemmas 1, 2, 4).** Lemma 1: term-by-term. Lemma 2/4
completeness: restrict an optimizer to its support I = supp λ*
(r = |I|); on R_{>0}^r the problem is smooth; if the constraint
gradients (1,…,1) and 2λ*_I are linearly dependent then λ*_I is
uniform (p₂ = 1/r, enumerated separately); otherwise Lagrange
multipliers exist and −ln λ − 1 = α + 2βλ forces ≤ 2 distinct
values. Enumerating two-level spectra over ranks and multiplicities
is therefore exhaustive for both the maximum and the minimum; the
extremal families are (a, rest-equal) for the max and (a×m₀, b×1),
m₀ = ⌊1/p₂⌋, for the min.

**B.2 (Lemma 3).** Size-biased variable P(X = λ_i) = λ_i: E X = p₂,
Var X = σ², δ = −E ln(1+u), u = X/p₂ − 1. Split at θp₂; on the large
part use −ln(1+u) ≤ −u + u²/(2θ²) (Taylor remainder with 1+ξ > θ);
on the small part bound Σλ ln(p₂/λ) at fixed mass w on ≤ d−1 atoms by
the equal split; bound w by Chebyshev, w ≤ σ²/((1−θ)²p₂²). Collecting
terms yields (5); minimizing the sum of the σ²/(2θ²p₂²) term and the
explicit inner maximum over θ gives (6) with limiting coefficient 1.

**B.3 (exact complete U-statistics and the UCB).** With ρ̄ the mean
marginal shadow on M (|M| qubits), Q = Σ_i ρ̂_i², and per-snapshot
invariants Tr ρ̂_i² = 5^{|M|}, Tr ρ̂_i³ = 7^{|M|}:

    p̂₂ = (N² Tr ρ̄² − N·5^{|M|}) / (N(N−1)),
    p̂₃ = (N³ Tr ρ̄³ − 3N Tr(Q ρ̄) + 2N·7^{|M|}) / (N(N−1)(N−2)),

verified against brute-force triple sums to 10⁻¹⁶. Exact influence
values: h₂,ᵢ = (N Tr(ρ̂_i ρ̄) − 5^{|M|})/(N−1), h₃,ᵢ =
Tr(ρ̂_i[(Nρ̄ − ρ̂_i)² − (Q − ρ̂_i²)])/((N−1)(N−2)); count-weighted
variances and covariance over the ≤ 6^{|M|} distinct outcome
patterns give Var(p̂₂) ≈ 4V₂/N, Var(p̂₃) ≈ 9V₃/N, Cov ≈ 6C₂₃/N, and
the delta-method variance of σ̂² = p̂₃ − p̂₂² + V̂ar(p̂₂). The UCB
uses two-sided (1 − α/10) normal intervals on the six active moment
bounds (Bonferroni) and evaluates (7) at the worst case over each
p₂-interval on a dense grid augmented by every flat point 1/r (where
the floors vanish — omitting these points would over-subtract and
break soundness; our first implementation did, and the bug was found
by exactly this argument).

**B.4 (variance model behind the extrapolations).** Certification
needs the added ceiling width 2 z ∂δ̄⁽³⁾/∂σ² · SE(σ̂²) ≲ τ −
U_population; with the §5 slope and z = 2.58 this is SE(σ̂²) ≲
5×10⁻⁴ at m = 4. The measured SE(σ̂²) at N = 10⁶ (m = 4, ghz-cl
states) is ≈ 4×10⁻⁴, matching the observed onset of certification;
SE scales as N^{−1/2} with a 2^{|M|}-type prefactor from the shadow
kernel variance, giving N*(witness, 80%) ≈ 2–3×10⁶ at m = 4 and
N* > 3^m for m ≤ 16 under the same model. These are model
extrapolations; the measured anchors are the table of §5.

## Appendix C: gate provenance

To prevent post-hoc threshold tuning in a negative benchmark, the
tolerance τ = 0.10 nats, the recall pool definition, noise grids,
recall targets, and the fallback branches were fixed before the runs
they judge (timestamped project log and seeds in the repository).
Recorded deviations: (1) the floor terms of Eq. (7) were added after
a two-seed pilot identified the entropy-gap cancellation as the
dominant noise cost, before the 10-seed validity map; (2) the
finite-shot estimator went through four recorded revisions (block
functional → subsampled U-statistics → covariance-aware box → exact
complete U-statistics), each replacing a diagnosed statistical flaw
(a non-vanishing Jensen bias of block-averaged nonlinear functionals;
a Monte-Carlo variance floor from subsample caps); revisions used
development data; a held-out coverage check on fresh seeds (350/350,
§5) closes the resulting in-sample concern; (3) the two pre-registered outcome branches for the
validity map did not anticipate a non-empty-but-smaller-than-targeted
certified region; the boundary framing used here was an explicit
decision recorded at that point; (4) a formally passing recall
clause based on a single separable state was identified and not
claimed; (5) a pre-run hypothesis that plug-in pipelines would
false-certify near threshold was not borne out (§5.1) and the
corresponding claim was downgraded under a rule recorded before that
run.

## Acknowledgments

The research design and the manuscript were stress-tested in multiple
rounds of adversarial review by large language models; findings —
including two estimator flaws, a proof regularity gap, the exact
U-statistic reduction of Appendix B.3, and several overclaim
corrections in this revision — were verified numerically before
adoption (transcripts and verification notes in the repository).
Part of the computations ran on Apple silicon; platform-specific
solver behavior is documented in the repository. All errors are ours.

## References

(Author lists, journal references and page numbers verified against
arXiv/APS on 2026-11-15; inline citations use the bracketed forms
above.)

- D. Petz, Commun. Math. Phys. 105, 123 (1986); Q. J. Math. 39, 97 (1988).
- P. Hayden, R. Jozsa, D. Petz, A. Winter, Commun. Math. Phys. 246, 359 (2004).
- O. Fawzi, R. Renner, Commun. Math. Phys. 340, 575 (2015). arXiv:1410.0664.
- M. Junge, R. Renner, D. Sutter, M. Wilde, A. Winter, Ann. Henri Poincaré 19, 2955 (2018). arXiv:1509.07127.
- M. Berta, K. Seshadreesan, M. Wilde, J. Math. Phys. 56, 022205 (2015). arXiv:1403.6102.
- N. Linden, M. Mosonyi, A. Winter, Proc. R. Soc. A 469, 20120737 (2013). arXiv:1212.0248.
- K. Seshadreesan, M. Wilde, Phys. Rev. A 92, 042321 (2015). arXiv:1410.1441.
- M. Berta, M. Tomamichel, IEEE Trans. Inf. Theory 62, 1758 (2016). arXiv:1502.07973.
- H.-Y. Huang, R. Kueng, J. Preskill, Nat. Phys. 16, 1050 (2020). arXiv:2002.08953.
- B. Vermersch, A. Rath, B. Sundar, C. Branciard, J. Preskill, A. Elben, PRX Quantum 5, 010352 (2024). arXiv:2304.12292.
- S. T. Flammia, Y.-K. Liu, Phys. Rev. Lett. 106, 230501 (2011). arXiv:1104.4695.
- Y.-A. Chen, C. Zhu, K. He, M. Jing, X. Wang, IEEE Trans. Inf. Theory 71, 5387 (2025). arXiv:2312.02031.
- P. Harremoës, F. Topsøe, IEEE Trans. Inf. Theory 47, 2944 (2001).
- D. W. Berry, B. C. Sanders, J. Phys. A 36, 12255 (2003). arXiv:quant-ph/0305059.
- J. Y. Lee, C.-M. Jian, C. Xu, PRX Quantum 4, 030317 (2023). arXiv:2301.05238.
- S. Sang, Y. Zou, T. H. Hsieh, Phys. Rev. X 14, 031044 (2024). arXiv:2310.08639.
- S. Sang, T. H. Hsieh, Phys. Rev. Lett. 134, 070403 (2025). arXiv:2404.07251.
- Y. Bao, R. Fan, A. Vishwanath, E. Altman, Phys. Rev. Lett. 136, 220402 (2026). arXiv:2301.05687.
- R. Fan, Y. Bao, E. Altman, A. Vishwanath, PRX Quantum 5, 020343 (2024). arXiv:2301.05689.
- R. Rubboli, A. Arqand, M. Wilde, arXiv:2609.34699.
- arXiv:2610.12929; arXiv:2610.11235 (stability of Markov length).
- Torvinen et al., arXiv:2605.06848.
- Petz-recovery experiments: W.-H. Png, V. Scarani, Phys. Rev. A 112, 022613 (2025), arXiv:2504.20399 (trapped ions); H. Li et al., arXiv:2606.12020 (photonic).
