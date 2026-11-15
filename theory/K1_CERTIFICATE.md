# K1 certificate theory (2026-11-13)

Conventions: all entropies in nats; H_α(ρ) = (1−α)^{-1} ln Tr ρ^α;
H1 = von Neumann; H2 = −ln p2, p2 = Tr ρ²; p3 = Tr ρ³. F_sq = ‖√ρ√σ‖₁².
Every formula below is re-derived here, not inherited (K0 lesson:
convention-pin everything; Carbon v5 formulas arrive in mixed
conventions).

## Lemma 1 (identity and one-sided bound)

For any tripartite ρ on S⊗F1⊗F2, with δ(M) := H1(M) − H2(M):

  I(S:F2|F1) − I'2 = δ(SF1) + δ(F1F2) − δ(F1) − δ(SF1F2),   (exact)
  0 ≤ δ(M)  always,
  ⇒ I(S:F2|F1) ≤ I'2 + δ(SF1) + δ(F1F2).                    (one-sided)

*Proof.* The identity is term-by-term: both I and I'2 are the same ±
combination of entropies of the four marginals, at α = 1 and α = 2
respectively. δ ≥ 0 is monotonicity of α ↦ H_α (Rényi). ∎

Chained with JRSWW (1509.07127, Thm 2.1, N = Tr_S; F_sq convention):

  −ln F_sq(ρ, R̄(ρ_{F1F2})) ≤ I(S:F2|F1) ≤ I'2 + δ(SF1) + δ(F1F2).

Only the two marginals SF1, F1F2 need certificates; spectral skew on F1
and SF1F2 only HELPS (their −δ ≤ 0 is discarded). Non-flat exact QMCs
(ρ⁺ anchor) make the inequality strict — false refusals, not unsoundness.

## Lemma 2 (k = 2 ceiling: exact, complete)

δ_max^(2)(p2, d) := max{H1(λ) : λ ∈ Δ_{d−1}, Σλᵢ² = p2} + ln p2.

(i) The maximizer has ≤ 2 distinct nonzero levels. *Proof.* Interior KKT
for max Σ−λlnλ s.t. Σλ = 1, Σλ² = p2: −lnλᵢ − 1 = α + 2βλᵢ for every
nonzero λᵢ, i.e. lnλ is an affine function of λ at all occupied levels;
ln λ (strictly concave) meets an affine function in ≤ 2 points. Zero
components reduce to the same problem at rank r < d. ∎

(ii) Therefore enumerating spectra (a×m, b×(r−m)) over 1 ≤ m < r ≤ d
(plus uniform rank-r when p2 = 1/r) is exhaustive — this is what
`delta_max2` does; machine-precision quadratic per (r, m).

(iii) Numerically the max coincides with the Harremoës–Topsøe family
x = (1 + √((d−1)(d·p2 − 1)))/d, rest equal (checked to 1e-15 on a probe
grid and pinned in tests as an independent implementation; we do not
rely on the literature claim — the enumeration is complete by (i)).

(iv) **k = 2 insufficiency (a result, not a footnote):**
δ_max^(2)(p2, d) = 0 iff p2 = 1 or p2 = 1/d. For intermediate-rank flat
states (p2 = 1/r, 1 < r < d) the ceiling is strictly positive — e.g.
δ_max^(2)(1/2, 4) = 0.2797 nats = 0.4035 bits: a GHZ-type marginal
cannot be certified flat from purity alone, at zero noise. This is what
forces p3 into the witness and the "low-order moments" (not
"second-moment") framing.

## Lemma 3 (k = 3 ceiling: sound closed-form bound)

Let σ² := p3 − p2² ≥ 0 (Cauchy–Schwarz on the size-biased variable;
σ² = 0 ⟺ λᵢ ∈ {0, p2} ∀i, i.e. the spectrum is flat at SOME rank).
For any θ ∈ (0, 1):

  δ(ρ) ≤ σ²/(2θ²p2²)
         + max_{0 ≤ w ≤ η(θ)} [ w·ln((d−1)p2 / w) − C(θ)·w ],
  η(θ) = min(1, σ²/((1−θ)²p2²)),  C(θ) = (1−θ) + (1−θ)²/(2θ²).

*Proof.* Size-biased distribution: P(X = λᵢ) = λᵢ; then E[X] = p2,
Var X = σ², and δ = E[ln(p2/X)] = −E[ln(1+u)], u := X/p2 − 1,
E[u] = 0, E[u²] = v := σ²/p2².

Split at θp2: S = {X ≤ θp2}, L = {X > θp2}; w := P(S) = Σ_{λ≤θp2} λ.
Since E[X] = p2 > θp2, L is nonempty, so S occupies ≤ d−1 atoms.

Large part: Taylor with remainder gives −ln(1+u) ≤ −u + u²/(2(1+ξ)²)
with ξ between 0 and u; on L, 1+ξ ≥ min(1, 1+u) > θ, hence
−ln(1+u) ≤ −u + u²/(2θ²). Taking expectations over L:
  E[−ln(1+u); L] ≤ −E[u; L] + E[u²; L]/(2θ²)
                 = E[u; S] + (v − E[u²; S])/(2θ²)
                 ≤ −(1−θ)w + (v − (1−θ)²w)/(2θ²),
using u ≤ θ−1 < 0 and u² ≥ (1−θ)² on S. This is v/(2θ²) − C(θ)w.

Small part: E[−ln(1+u); S] = Σ_S λ ln(p2/λ). At fixed mass Σ_S λ = w on
≤ d−1 atoms the concave sum is maximized by the equal split λ = w/(d−1):
value w·ln((d−1)p2/w) (the cap λ ≤ θp2 can only reduce it).

Mass bound: w = P(p2 − X ≥ (1−θ)p2) ≤ Var X/((1−θ)²p2²) (Chebyshev),
and w ≤ 1. Maximizing over admissible w and minimizing over θ gives the
statement. ∎

Properties: bound = 0 on every flat spectrum of any rank (σ² = 0 ⇒ both
terms vanish); continuous in (p2, p3); any single θ is already sound —
the θ-grid in `delta_max3_bound` affects tightness only. The inner max
is explicit: g(w) = w ln(A/w) − Cw is concave with g(w*) = w* at
w* = A·e^{−1−C}, so max = g(min(w*, η)).

## Lemma 4 (≤ 3 levels at fixed p2 AND p3)

Interior KKT for max H1 s.t. (Σλ, Σλ², Σλ³) fixed:
−lnλ − 1 = α + 2βλ + 3γλ² at occupied levels: lnλ equals a quadratic in
λ in ≤ 3 points (strict concavity of ln vs parabola). So the exact
δ_max^(3) is attained on ≤ 3-level spectra. `delta_max3_exact`
enumerates them with multi-start root finding — BEST-EFFORT ONLY: a
missed root UNDERESTIMATES the ceiling, so the exact search is never
the witness; it measures the slack of Lemma 3. The sound witness is
always the Lemma-3 bound.

## Accessibility (honest costs; estimator design is K3 scope)

p2 and p3 of SF1 and F1F2 are functionals of the same single-qubit
Pauli bitstrings (2nd/3rd-order U-statistics across snapshots); the p3
estimator's variance grows as O(4^|M|) vs O(2^|M|) for p2 — the witness
is "cheap" only while |SF1|, |F1F2| stay moderate. Finite-shot
one-sided confidence bounds (δ = 0.05) replace point estimates in K3.
[to-verify before preprint: the standard reference set for multi-copy
moment estimators from randomized measurements — candidate sources
van Enk–Beenakker 2012, Elben et al. toolbox — verify exact claims
before citing; no citation is load-bearing here.]

## What this closes and what it leaves open (O27 residue)

Closed: soundness of U^(k) as an upper bound on −ln F_sq(ρ, R̄) for
k ∈ {2, 3}; exactness of the k=2 ceiling; zero-on-flat and continuity
of the k=3 bound; the k=2 insufficiency result.
Open (owed to K2): recall — how much of the recoverable world the
certified region covers under noise; the measured slack of Lemma 3 vs
the exact 3-level ceiling; behaviour on the HJPW false-refusal family.

## Lemma 5 (δ_min^(2); added 2026-11-13 post-K2-smoke, pre-full-map)

δ_min^(2)(p2, d) := min{H1(λ) : λ ∈ Δ_{d−1}, Σλᵢ² = p2} + ln p2 ≥ 0 is
computed EXACTLY by the same two-level enumeration as Lemma 2 (the
≤2-level KKT argument constrains every interior extremum, minima
included; boundary recurses). Regularity (v6 patch, 2026-11-14):
restricting an optimizer to its support I = supp(λ*) (r = |I|) gives an
interior extremum on R_{>0}^r; either λ*_I is uniform (1 level,
p2 = 1/r — enumerated separately), or the constraint gradients
∇Σμᵢ and ∇Σμᵢ² are linearly independent (they are parallel only at the
uniform point), so Lagrange multipliers exist and the ≤2-root argument
applies. Second-order refinement (verified to machine precision against
the enumeration): the minimum is attained at (a×m, b×1) with
m = ⌊1/p2⌋, a = (m+√(m((m+1)p2−1)))/(m(m+1)), b = (1−√(m((m+1)p2−1)))/(m+1)
— attribution: Harremoës–Topsøe, IEEE TIT 47, 2944 (2001), Thm II.5;
Berry–Sanders, J. Phys. A 36, 12255 (2003), Sec. 2 [verify theorem/page
numbers at draft time]. Our enumeration remains the implementation (it
is complete regardless). Corollary (improved witness):

  I(S:F2|F1) ≤ U^(3)± := I'2 + δ_max-bound(SF1) + δ_max-bound(F1F2)
                         − δ_min^(2)(F1) − δ_min^(2)(SF1F2),

directly from the Lemma-1 identity; δ_min at fixed (p2,p3) ≥ δ_min at
fixed p2, so the k=2 minimum is a valid certified lower bound.

Provenance note (honesty box material): the K2 smoke run (2 seeds)
diagnosed that the one-sided witness's dominant noise cost is the
DISCARDED cancellation −δ(F1)−δ(SF1F2) (measured on noisy classical
GHZ: δ(SF1)+δ(F1F2) = 0.0746 vs discarded 0.0713 at CMI = 0.0033,
p = 5e-3). U^(3)± was registered after that diagnosis and before the
full 10-seed map; the ratified K2 gate is still judged on U^(2)/U^(3),
with U^(3)± reported alongside as the registered improvement.
