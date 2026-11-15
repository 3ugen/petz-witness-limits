# One-sided low-moment witnesses of twirled-Petz recoverability: noise validity boundaries and shot-budget limits against Pauli state reconstruction

**Draft v0.1 (2026-11-14). Target: PRA / Quantum (methods and limits).**

## Abstract

Deciding whether a tripartite state ρ_{SF₁F₂} can be approximately
recovered from its marginal ρ_{F₁F₂} by a channel acting on F₁ alone is
controlled by the conditional mutual information I(S:F₂|F₁) through the
Fawzi–Renner/JRSWW bound −ln F²(ρ, R̄(ρ_{F₁F₂})) ≤ I, with R̄ the
universal twirled Petz map. We ask what the cheapest measurable
quantities — second and third Pauli moments of marginals — can certify
about this recoverability, and at what price. We construct a one-sided
witness U = I′₂ + δ̄(SF₁) + δ̄(F₁F₂) − δ̲(F₁) − δ̲(SF₁F₂), where I′₂ is
the (provably CMI-unlike) Rényi-2 entropy combination and δ̄, δ̲ are
exact or closed-form ceilings/floors on the entropy gap H₁ − H₂ at
fixed low moments; soundness I ≤ U is an exact identity plus moment
extremization, never an approximation. We prove that purities alone
cannot certify any intermediate-rank flat state (an irreducible
Θ(log d) ceiling), map the witness's validity on 9,680 simulated
state–noise points (zero soundness violations; certified-region noise
boundary p_c(τ = 0.1 nats) ≈ 3×10⁻³ per-qubit depolarizing, extending
to ≈10⁻² under dephasing-biased noise, with its mechanism — the
forfeited cancellation of marginal entropy gaps — measured directly),
and benchmark finite-shot certification head-to-head against projected
MLE reconstruction on identical random-Pauli data. The budget verdict
is negative and we report it as the main empirical result: at m = 4
qubits the witness certifies only from ~10⁶ snapshots while
reconstruction reaches the same decisions from ~3×10⁵, and neither
method certifies in the genuinely sub-tomographic regime
N_settings < 3^m for m ≥ 4 at practical budgets. What survives is the
guarantee structure: the witness's refusals carry an order-of-magnitude
larger safety margin on near-threshold negatives (minimum UCB−τ of
+0.45 vs +0.05 nats for plug-in SDP fidelity-of-recovery), and its
soundness is a theorem rather than an empirical observation. We state
precisely when "just reconstruct" is the right advice, and what a
moment-based certificate is still uniquely good for.

## 1. Introduction

Approximate quantum Markov chains sit at the junction of several
currently active lines: Petz recovery maps as near-optimal decoders,
mixed-state phases diagnosed by information quantities on decohered
stabilizer states, and the operational reading of quantum Darwinism as
recoverability of a system from environment fragments [Torvinen et al.,
arXiv:2605.06848]. The quantitative anchor is the Fawzi–Renner theorem
and its universal-recovery refinement [Fawzi–Renner, CMP 340, 575
(2015), arXiv:1410.0664; Junge–Renner–Sutter–Wilde–Winter, Ann. Henri
Poincaré 19, 2955 (2018), arXiv:1509.07127]: for the twirled Petz map
R̄ built from ρ_{SF₁} and ρ_{F₁},

    −ln F²(ρ_{SF₁F₂}, R̄_{F₁→SF₁}(ρ_{F₁F₂})) ≤ I(S:F₂|F₁).            (1)

(Throughout, F² is squared/Uhlmann fidelity and all information
quantities are in nats; Fawzi–Renner's F is root fidelity, which
accounts for the absent factor 2.)

Certifying recoverability through (1) requires an upper bound on the
CMI, and the CMI is expensive: its von Neumann entropies are not
functions of any finite set of polynomial moments. Experiments that
implement Petz recovery to date verify it by full tomography
[arXiv:2606.12020; arXiv:2504.20399]. At the other extreme, the
cheapest measurable quantities on Pauli data are low moments
Tr ρ_M^k of marginals, and the Rényi-2 combination

    I′₂ = S₂(SF₁) + S₂(F₁F₂) − S₂(F₁) − S₂(SF₁F₂)                    (2)

is emphatically NOT a conditional mutual information: Berta,
Seshadreesan and Wilde explicitly disqualify such linear combinations
as Rényi generalizations of CMI [JMP 56, 022205 (2015),
arXiv:1403.6102], and Linden–Mosonyi–Winter prove there are no
universal linear inequalities among Rényi marginal entropies at all
[Proc. R. Soc. A 469, 20120737 (2013), arXiv:1212.0248]. We exhibit
the pathology concretely: two classical Markov chains (perfectly
Petz-recoverable, CMI = 0) with I′₂ = ln(8/9) < 0 and
I′₂ = ln(10/9) > 0 respectively (§3.1). In decohered stabilizer
phases the same gap has physical teeth: Rényi-2 diagnostics see their
own transition (p_c^(2) ≈ 0.178 for the toric code under bit flips)
while true decodability is lost already at p_c ≈ 0.109 [Bao–Fan–
Vishwanath–Altman, arXiv:2301.05687; Fan–Bao–Altman–Vishwanath,
arXiv:2301.05689] — a raw I′₂ ≈ 0 in the window between the two
thresholds falsely signals decodability. Any honest moment-based
recoverability statement must therefore be one-sided and must carry
its own correction terms.

Recent Rényi-CMI proposals — including the log-Euclidean Rényi CMI of
Rubboli, Arqand and Wilde [arXiv:2609.34699] — preserve recovery
bounds at the density-matrix level via matrix powers and logarithms
that require full spectral access and are not functions of finite
polynomial moments. Conversely, virtual quantum Markov chain SDPs
[Chen et al., IEEE TIT 71, 5387 (2025), arXiv:2312.02031] certify
HPTP ("virtual") rather than physical CPTP recovery, and do so from a
known classical description of the state, not from data; their
Example 4 (equal-CMI states with different virtual recoverability)
shows moment-based quantities cannot reach the HPTP level at all. The
stability of CMI decay itself under local noise is a live question
[arXiv:2610.12929; arXiv:2610.11235]. Here we ask the complementary
measurement-level question: **when can second and third Pauli moments
alone provide a sound one-sided certificate of physical twirled-Petz
recoverability, which noise mechanism sets the certificate's recall
boundary even when the true CMI remains near zero, and how do its
finite-shot budgets compare, on identical data, against simply
reconstructing the state?**

Our answers, in brief:

1. **A sound certificate exists and is closed-form** (§3): the exact
   identity I − I′₂ = δ(SF₁) + δ(F₁F₂) − δ(F₁) − δ(SF₁F₂), with
   δ(M) = H₁(M) − H₂(M) ≥ 0, turns into I ≤ U after replacing the
   positive pair by moment ceilings and the negative pair by moment
   floors — exact at k = 2 moments (complete two-level extremization,
   recovering Harremoës–Topsøe / Berry–Sanders single-marginal forms)
   and a continuous closed-form bound at k = 3 via the size-biased
   spectral variance σ² = p₃ − p₂², which vanishes iff the spectrum is
   flat at any rank.
2. **Purities alone provably do not suffice** (§3.2): the k = 2
   ceiling vanishes only at p₂ ∈ {1, 1/d}; a GHZ-type rank-2 marginal
   carries an irreducible 0.28-nat ceiling at zero noise. "Second
   moments" is not a viable pitch; "second and third" is.
3. **The validity map** (§4): on 9,680 state–noise points across eight
   seeded ensembles the witness commits zero soundness violations; its
   certified region is non-empty and perfectly reproducible, with
   noise boundary p_c(τ = 0.10) = 3×10⁻³ per-qubit depolarizing on
   stabilizer-Markov states — set not by the bound's slack but by the
   forfeited cancellation between marginal entropy gaps, which we
   measure directly. The boundary is channel-dependent (≈10⁻² under
   10:1 dephasing-biased noise) and scales roughly linearly in τ.
4. **The budget verdict is negative** (§5): with exact complete
   U-statistic estimators (no subsampling), certification at m = 4
   needs ~10⁶ random-Pauli snapshots, while projected-MLE
   reconstruction of the full 16-dimensional state reaches 80% recall
   at the same τ from 10⁶ and 50% already at 3×10⁵ — reconstruction
   dominates by a factor ~2–3 wherever either method works, and
   neither certifies in the sub-tomographic regime N_settings < 3^m.
   Below m ≈ 10, "just reconstruct" is the right advice, and we say so.
5. **What survives** (§6): the witness's refusals are a theorem, not an
   observation. On near-threshold non-recoverable states (exact
   infidelity 0.13–0.18 at τ = 0.10), 100-repetition runs show zero
   false certifications for all pipelines, but the witness's minimum
   refusal margin is +0.45 nats versus +0.05 for plug-in SDP
   fidelity-of-recovery — an order of magnitude. Where a wrong
   "recoverable" verdict is expensive, the guarantee structure, not
   the shot count, is the deliverable.

Everything in this paper was executed under pre-registered KPI gates
with recorded NO-GO branches; design changes made after seeing data
are flagged as such (Appendix C, "honesty box").

## 2. Setup and conventions

Tripartition S ⊗ F₁ ⊗ F₂ with dimensions (d_S, d₁, d₂); all numerics
use qubit registers with contiguous blocks. Fidelity is squared/
Uhlmann, F²(ρ,σ) = ‖√ρ√σ‖₁²; logs are natural. H_α denotes Rényi
entropy, H₁ von Neumann, H₂ = −ln p₂ with p₂ = Tr ρ², p₃ = Tr ρ³.
The twirled Petz map is R̄(·) = ∫ β₀(t) R^{t/2}(·) dt with
β₀(t) = (π/2)(cosh πt + 1)^{-1} and rotated maps R^{t/2} as in JRSWW;
in the eigenbases of ρ_{SF₁} and ρ_{F₁} its matrix elements carry the
closed-form kernel ω/sinh ω at half-log phase ω = ½ ln(λ_a μ_y / λ_b
μ_x) (we verified this kernel against direct quadrature of β₀; a
convention slip between half-log and full-log phases is easy to make
and was caught in our own pipeline by the two-path cross-check). Maps
are completed to CPTP outside supp ρ_{F₁}; all oracle implementations
are duplicated through independent code paths agreeing to 10⁻¹⁰.

## 3. The certificate

### 3.1 Identity and one-sided bound

**Lemma 1.** For any ρ on S ⊗ F₁ ⊗ F₂, with δ(M) := H₁(M) − H₂(M) ≥ 0:

    I(S:F₂|F₁) − I′₂ = δ(SF₁) + δ(F₁F₂) − δ(F₁) − δ(SF₁F₂),          (3)

exactly; hence I ≤ I′₂ + δ(SF₁) + δ(F₁F₂), and combining with (1),

    −ln F²(ρ, R̄(ρ_{F₁F₂})) ≤ I ≤ I′₂ + δ(SF₁) + δ(F₁F₂).            (4)

The proof is term-by-term (both sides of (3) are the same ±combination
of marginal entropies at α = 1 and α = 2) plus monotonicity of
α ↦ H_α. Only SF₁ and F₁F₂ need upper certificates; the discarded
−δ(F₁) − δ(SF₁F₂) ≤ 0 can only strengthen (4) — or be partially
reclaimed, which is Lemma 5.

The classical anchors ρ⁻ and ρ⁺ (three-bit Markov chains; explicit
atoms in Appendix A) have CMI = 0, exact Petz recovery fidelity 1, and
I′₂ = ln(8/9) and ln(10/9): sign-indefiniteness of I′₂ on perfectly
recoverable states is not a pathology of exotic states.

### 3.2 Moment ceilings and floors

δ(M) is not measurable, but its extreme values at fixed low moments
are computable:

**Lemma 2 (k = 2 ceiling, exact and complete).** The maximizer of H₁
at fixed (p₂, d) has at most two distinct nonzero eigenvalues (interior
Lagrange stationarity forces ln λ affine in λ; zero eigenvalues recurse
to lower rank; at the only point where the constraint gradients are
parallel the spectrum is uniform and enumerated separately).
Enumerating two-level spectra is therefore exhaustive; the maximum is
the Harremoës–Topsøe one-high-rest-flat family [IEEE TIT 47, 2944
(2001); Berry–Sanders, J. Phys. A 36, 12255 (2003)]. We write
δ̄⁽²⁾(p₂, d) for the resulting exact ceiling.

**Proposition (k = 2 insufficiency).** δ̄⁽²⁾(p₂, d) = 0 iff p₂ = 1 or
p₂ = 1/d. For every intermediate-rank flat state (p₂ = 1/r, 1 < r < d)
the ceiling is strictly positive: δ̄⁽²⁾(1/2, 4) = 0.2797 nats (0.4035
bits). A GHZ-type marginal cannot be certified flat from purity alone,
at zero noise. This is why the witness needs p₃.

**Lemma 3 (k = 3 ceiling, closed-form and sound).** With
σ² := p₃ − p₂² ≥ 0 (the variance of the size-biased spectral variable;
σ² = 0 iff the spectrum is flat at some rank), for every θ ∈ (0,1)

    δ ≤ σ²/(2θ²p₂²) + max_{0≤w≤η} [ w ln((d−1)p₂/w) − C(θ) w ],      (5)

with η = min(1, σ²/((1−θ)²p₂²)) and C(θ) = (1−θ) + (1−θ)²/(2θ²); the
derivation (size-biased split, Taylor remainder, Chebyshev) is in
Appendix B. The bound is zero on every flat spectrum of any rank and
continuous in (p₂, p₃); minimizing over θ, its small-σ² asymptotics is

    δ̄⁽³⁾(σ²; p₂, d) ≈ 1.4 · (σ²/p₂²) · ln((d−1)p₂³/σ²)              (6)

(the prefactor is our measured envelope over σ² ∈ [10⁻⁵, 10⁻²]; the
slope is log-divergent as σ² → 0, a fact that controls the shot
budgets of §5).

**Lemma 5 (k = 2 floor, exact).** The same two-level enumeration
computes δ̲⁽²⁾(p₂, d) = min δ at fixed p₂ (attained at m = ⌊1/p₂⌋
copies of one level plus a single smaller one). Since the minimum at
fixed (p₂,p₃) is no smaller, subtracting δ̲⁽²⁾ of the discarded
marginals preserves (4):

    I ≤ U := I′₂ + δ̄⁽³⁾(SF₁) + δ̄⁽³⁾(F₁F₂) − δ̲⁽²⁾(F₁) − δ̲⁽²⁾(SF₁F₂). (7)

U is the witness used throughout; "certify" means U ≤ τ with the
pre-registered tolerance τ = 0.10 nats (F² ≥ 0.905 for R̄).

The Bao-window observation of §1 now reads constructively: in the
(0.109, 0.178) regime of decohered stabilizer codes where raw I′₂ ≈ 0
is unsound, the ceilings in (7) are large — U refuses. The correction
terms are not decoration; they are what makes a Rényi-2-based quantity
usable at all near such windows.

## 4. Validity map

We evaluated (7) with exact moments against exact oracles on 9,680
points: eight seeded ensembles (stabilizer-Markov states with
entangling perturbations; bounded-rank random; central-spin dephasing;
classical branching with pointer bias; 1D Gibbs chains; non-flat exact
Markov chains of Hayden–Jozsa–Petz–Winter form [CMP 246, 359 (2004)];
W–GHZ interpolations; a repetition-code family) × a per-qubit
depolarizing grid × 10 seeds, on 4–5-qubit tripartitions. Headlines:

- **Soundness: 0 violations in 9,680 points** (as the theorem demands;
  the map doubles as a regression of the entire oracle stack).
- **Recall on recoverable, non-trivially-correlated states** (−ln F² ≤
  0.02 and I(S:F₁) ≥ ½ ln 2): on the stabilizer-Markov ensemble, 69%
  at p = 0, 59% at 10⁻³, 36% at 3×10⁻³, 8% at 5×10⁻³, 0 beyond. The
  purity-only witness U⁽²⁾ certifies **nothing anywhere** — the §3.2
  insufficiency is a map-wide empirical fact, not only a theorem.
- **p_c(τ = 0.10) = 3×10⁻³** per-qubit depolarizing, with zero spread
  across 10 ensemble seeds; approximately linear in τ (1.1×10⁻³ at
  τ = 0.05, 10⁻² at τ = 0.25); channel-dependent: under 10:1
  dephasing-biased Pauli noise the certified region survives to
  ≈10⁻² (recall 35% at 5×10⁻³, 21% at 10⁻²), while amplitude damping
  tracks depolarizing. All p_c statements in this paper carry their
  (τ, m, channel) parameters.
- **The mechanism is measured, not conjectured**: on noisy classical
  GHZ at p = 5×10⁻³ the true CMI is 0.0033 and the recovery
  infidelity 0.0024 — deep in the recoverable regime — yet
  δ(SF₁) + δ(F₁F₂) = 0.0746 while the discarded pair δ(F₁) + δ(SF₁F₂)
  = 0.0713. The one-sided witness pays the full noise-induced
  flatness defect of two marginals and forfeits a near-perfect
  cancellation; the δ̲-floors of Lemma 5 reclaim roughly a third of
  it (p_c: 2×10⁻³ → 3×10⁻³). This cost is intrinsic to any one-sided
  moment certificate built on (3), and it bites while the true CMI is
  still negligible — a sharper statement than bound slack.
- **Structural false refusals**: Gibbs chains and non-flat HJPW exact
  Markov chains are recoverable but never certified (recall 0 even at
  p = 0) — soundness is preserved, completeness is the designed
  sacrifice. The witness detects "flat-spectrum Markovianity", and we
  claim nothing broader.

## 5. Finite-shot budgets: witness vs reconstruction

Estimation uses random-Pauli classical shadows [Huang–Kueng–Preskill,
Nat. Phys. 16, 1050 (2020)]: every snapshot is one shot in an i.i.d.
uniform product Pauli basis. Both pipelines consume the SAME
bitstrings:

- **Witness pipeline**: exact complete U-statistics for p₂ and p₃ of
  the four marginals via a Newton–Girard reduction on the accumulated
  marginal shadow and its square (closed forms in Appendix B.3; cost
  O(min(N, 6^|M|)·d_M²), no subsampling), exact first-order influence
  functions for variances and the p̂₂–p̂₃ covariance, the bias
  correction σ̂² = p̂₃ − p̂₂² + V̂ar(p̂₂), and a Bonferroni moment-box
  propagated through the monotone structure of (7) into a one-sided
  95% upper confidence bound Û. (Degenerate marginals — e.g. a
  maximally mixed F₁, where the first-order projection vanishes —
  make the normal quantile approximate; the UCB is asymptotic
  (1−δ) on non-degenerate states and empirically covered the exact
  value in 100% of our runs.)
- **Reconstruction pipeline**: projected MLE (PSD, unit-trace) of the
  full ρ_{SF₁F₂} from the same shadows, plug-in CMI, delete-block
  jackknife 95% UCB. At m = 4 the budgets below put hundreds of shots
  in each of the 3⁴ = 81 bases — this is full-basis Pauli state
  reconstruction, and we label it as such.

**Result (m = 4, dims 2×4×2; six states, 6 repetitions/point):**

| N snapshots | witness recall | reconstruction recall |
|---|---|---|
| 10⁴ | 0 | 0 |
| 3×10⁴ | 0 | 12.5% |
| 10⁵ | 0 | 37.5% |
| 3×10⁵ | 0 | 50% |
| 10⁶ | 29.2% | 100% |

Reconstruction reaches 80% recall at N* = 10⁶; the witness crosses
zero recall only at 10⁶ (extrapolated N* ≈ 2–3×10⁶). UCB coverage of
the exact value: 100% (witness), 98.3% (reconstruction); zero false
certifications of negatives for both. At m = 6, neither method
certifies anything up to N = 3×10⁵, including the genuinely
sub-tomographic point N = 300 < 3⁶; a variance-scaling estimate from
(6) puts the witness's N* above 3^m for all m ≤ 16 at τ = 0.10.

The driver is (6): the ceiling's log-divergent sensitivity to σ̂² means
certification at headroom τ needs SE(σ̂²) of order τ/(70·…) — about
6×10⁻⁴ at m = 4 — while reconstruction of a d = 16 state is cheap in
absolute terms. **For marginal sizes where either method is practical,
reconstruction wins on shots by ×2–3; below m ≈ 10, "just
reconstruct" is the right engineering advice.** We looked for a
regime where the witness is cheaper and did not find one; this
negative is a primary result, pre-registered as such.

### 5.1 Near-threshold negatives: what the guarantee buys

On the mixture family (1−λ)|GHZ₄⟩⟨GHZ₄| + λ ρ_GHZ-cl tuned to exact
twirled-Petz infidelity 0.13 and 0.18 (true negatives just above
τ = 0.10), 100 repetitions per budget N ∈ {3×10³, 10⁴, 3×10⁴}:

- **False certifications: zero, for all three pipelines** (witness
  0/600; MLE-CMI 0/600; plug-in SDP fidelity-of-recovery
  [Seshadreesan–Wilde, PRA 92, 042321 (2015); SDP form after
  Berta–Tomamichel, IEEE TIT 62, 1758 (2016)] 0/60). The anticipated
  decision-level unsoundness of plug-ins did not materialize at these
  budgets, and we report that as measured.
- **Refusal margins differ by an order of magnitude**: minimum
  UCB − τ across all runs: witness +0.45 / +0.54 nats (two targets),
  MLE-CMI +0.19 / +0.21, FoR-SDP +0.05 / +0.08. Point-estimate probes
  on a far-threshold negative show the known downward plug-in bias
  (FoR-SDP estimates infidelity 0.43 where the truth is 0.693 at
  N = 3×10⁴).

The honest summary: plug-in pipelines were empirically safe here but
run close to the edge with no guarantee; the witness refuses with a
wide, theorem-backed margin. Where a false "recoverable" is expensive
— decoder qualification, erasure-resilience claims — that asymmetry,
not the shot count, is what a moment certificate is for.

## 6. Discussion

**When to use what.** If you can afford N ≳ 4^m informationally
complete snapshots and tolerate guarantee-free error bars:
reconstruct, plug in, jackknife — it is cheaper and it was never
caught false-certifying in our runs. Use the moment witness when (i)
the decision is one-sided and safety-critical, (ii) the states of
interest are near the flat-spectrum manifold (stabilizer-grade
preparations at per-qubit error ≲ 3×10⁻³, or ≲10⁻² under
dephasing-dominant noise), and (iii) a closed-form, auditable
certificate is worth more than shot economy. The certificate's
by-products — the exact k = 2 ceilings/floors and the σ²-flatness
characterization — are independently useful wherever "is this marginal
flat?" is the real question.

**Relation to nearby work.** Our object is physical CPTP (twirled-
Petz) recoverability from data; virtual-recovery SDPs [2312.02031]
answer a different (HPTP, known-state) question that moments provably
cannot see. Log-Euclidean and other proper Rényi CMIs [2609.34699;
1403.6102] carry recovery bounds but are not moment-accessible; our
U is deliberately NOT a Rényi CMI — it is an engineered envelope
around the non-CMI quantity I′₂. The decohered-stabilizer literature
[2301.05238; 2310.08639; 2404.07251; 2301.05687; 2301.05689] supplies
both the states where the witness is at its best and the cautionary
two-thresholds picture that motivates its correction terms. For
quantum Darwinism, the witness gives a cheap sufficient check of the
Petz-plateau reading [2605.06848] on fragment pairs in its validity
region.

**Open directions.** (i) Many-pairs amortization: across families of
overlapping (F₁, F₂) tripartitions the same snapshots feed every
marginal, and the witness's relative economics improve; we did not
benchmark this axis. (ii) Quantum marginal constraints: the moment
boxes treat (p₂, p₃) of different marginals as independent, which is
loose. (iii) Two-copy access: σ² is a two-copy observable; swap-test
data would bypass the 2^m variance wall and reopen the budget
question. (iv) The instability of CMI decay under local channels
[2610.12929] suggests mapping p_c for witness classes beyond
flat-spectrum.

## Appendix A: anchors and ensembles

ρ⁻ = ½|00⟩⟨00|_{SF₁} ⊗ I/2 + ½ (I/2)_S ⊗ |10⟩⟨10|_{F₁F₂};
ρ⁺ = ½|000⟩⟨000| + ½ (I/2)_S ⊗ |1⟩⟨1|_{F₁} ⊗ (I/2)_{F₂}. Both are
classical Markov chains (CMI = 0, Petz fidelity 1) with
I′₂ = ln(8/9) and ln(10/9). Ensemble definitions, seeds, and the full
noise grids are in the repository (ensembles.py; every figure and
table regenerates from seeded scripts).

## Appendix B: proofs and estimator formulas

B.1 Lemmas 1–2, 5 as in §3 with the regularity discussion (LICQ fails
only at uniform spectra, which are enumerated separately). B.2 Proof
of (5): under the size-biased distribution P(X = λ_i) = λ_i one has
E X = p₂, Var X = σ², δ = −E ln(1+u) with u = X/p₂ − 1; split at θp₂,
bound the large part by second-order Taylor with remainder ≥ θ, the
small part by concavity at fixed mass on ≤ d−1 atoms, and the mass by
Chebyshev. B.3 Exact complete U-statistics: with ρ̄ the mean marginal
shadow, Q = Σ_i ρ̂_i², and per-snapshot invariants Tr ρ̂_i² = 5^k,
Tr ρ̂_i³ = 7^k,

    p̂₂ = (N² Tr ρ̄² − N·5^k) / (N(N−1)),
    p̂₃ = (N³ Tr ρ̄³ − 3N Tr(Q ρ̄) + 2N·7^k) / (N(N−1)(N−2)),

with exact influence values h₂,ᵢ, h₃,ᵢ in closed form (repository,
sampling.py); verified against brute-force triple sums to 10⁻¹⁶.

## Appendix C: honesty box (methods provenance)

Pre-registration: all gates, thresholds (τ = 0.10, recall targets,
noise grids) and NO-GO branches were recorded before the runs they
judge; the project log is in the repository. Deviations, flagged when
made: (1) the δ̲-floor upgrade (Lemma 5 / Eq. 7) was registered after a
2-seed smoke run diagnosed the cancellation mechanism and before the
full 10-seed map; (2) the finite-shot estimator went through four
versions (block-functional → incomplete-U box → covariance-aware →
exact Newton–Girard), each replaced after its specific statistical
flaw was diagnosed (Jensen bias of block-averaging; a Monte-Carlo
variance floor from subsample caps); all versions and their record
runs are preserved; (3) the K2 gate's two pre-registered branches
(GO/NO-GO) did not anticipate a non-empty-but-smaller-than-targeted
certified region; the "boundary" reading used here was an explicit
owner decision, recorded; (4) a formally passing GO clause via a
single separable state was identified and NOT claimed; (5) the
anticipated decision-level unsoundness of plug-in baselines did not
appear in the near-threshold run, and the corresponding claim was
downgraded per a rule recorded before that run. Numerical cross-checks
caught, during development: a twirl-kernel convention slip (half-log
vs full-log phase), a partial-trace axis-pairing bug, a MINPACK
false-convergence mode (ier = 1 at a non-root under singular
Jacobians), and a CLARABEL SolverError class on rank-deficient
reconstructions (macOS arm64) with an SCS fallback now in code.

## Acknowledgments

The research design was stress-tested through multiple rounds of
adversarial review by large language models (Google Gemini "Carbon"
series), whose findings — including two estimator flaws, one proof
regularity gap, and the exact-U-statistic reduction of Appendix B.3 —
were independently verified numerically before adoption; reviews and
verification notes are preserved in the repository. Part of the
computations ran on Apple silicon; platform-specific solver behavior
is documented in the honesty box. All errors are ours.

## References

(Repository-verified set; page/theorem numbers to be re-checked at
submission.)

- D. Petz, Commun. Math. Phys. 105, 123 (1986); Q. J. Math. 39, 97 (1988).
- P. Hayden, R. Jozsa, D. Petz, A. Winter, Commun. Math. Phys. 246, 359 (2004).
- O. Fawzi, R. Renner, Commun. Math. Phys. 340, 575 (2015). arXiv:1410.0664.
- M. Junge, R. Renner, D. Sutter, M. Wilde, A. Winter, Ann. Henri Poincaré 19, 2955 (2018). arXiv:1509.07127.
- D. Sutter, R. Renner, Lett. Math. Phys. 108, 2007 (2018). arXiv:1705.06749.
- D. Sutter, Approximate Quantum Markov Chains, Springer (2018). arXiv:1907.09421.
- M. Berta, K. Seshadreesan, M. Wilde, J. Math. Phys. 56, 022205 (2015). arXiv:1403.6102.
- N. Linden, M. Mosonyi, A. Winter, Proc. R. Soc. A 469, 20120737 (2013). arXiv:1212.0248.
- K. Seshadreesan, M. Wilde, Phys. Rev. A 92, 042321 (2015). arXiv:1410.1441.
- M. Berta, M. Tomamichel, IEEE Trans. Inf. Theory 62, 1758 (2016). arXiv:1502.07973.
- H.-Y. Huang, R. Kueng, J. Preskill, Nat. Phys. 16, 1050 (2020). arXiv:2002.08953.
- B. Vermersch et al., PRX Quantum 5, 010352 (2024). arXiv:2304.12292.
- Y.-A. Chen, C. Zhu, K. He, M. Jing, X. Wang, IEEE Trans. Inf. Theory 71, 5387 (2025). arXiv:2312.02031.
- P. Harremoës, F. Topsøe, IEEE Trans. Inf. Theory 47, 2944 (2001).
- D. W. Berry, B. C. Sanders, J. Phys. A 36, 12255 (2003). arXiv:quant-ph/0305059.
- J. Y. Lee, C.-M. Jian, C. Xu, PRX Quantum 4, 030317 (2023). arXiv:2301.05238.
- S. Sang, Y. Zou, T. H. Hsieh, Phys. Rev. X 14, 031044 (2024). arXiv:2310.08639.
- S. Sang, T. H. Hsieh, Phys. Rev. Lett. 134, 070403 (2025). arXiv:2404.07251.
- Y. Bao, R. Fan, A. Vishwanath, E. Altman, arXiv:2301.05687.
- R. Fan, Y. Bao, E. Altman, A. Vishwanath, arXiv:2301.05689.
- R. Rubboli, A. Arqand, M. Wilde, arXiv:2609.34699.
- arXiv:2610.12929 (instability of quantum Markov length); arXiv:2610.11235.
- Torvinen et al., arXiv:2605.06848.
- Petz-recovery experiments: arXiv:2504.20399 (trapped ions); arXiv:2606.12020 (photonic); arXiv:2510.26895 (tabletop reversibility conditions).
- E. H. Lieb, M. B. Ruskai, J. Math. Phys. 14, 1938 (1973).
