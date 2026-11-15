"""K1 certificate: moment ceilings on the Rényi gap δ(M) = H1(M) − H2(M).

Theory (proofs in theory/K1_CERTIFICATE.md, convention-pinned per the K0
lesson; everything in NATS):

- Lemma 1 (identity + one-sided bound): I − I'2 = δ(SF1)+δ(F1F2)
  −δ(F1)−δ(SF1F2) exactly, with δ(M) ≥ 0 by monotonicity of α ↦ H_α, so
  I ≤ I'2 + δ(SF1) + δ(F1F2) unconditionally.
- Lemma 2 (k = 2 ceiling, EXACT and complete): at fixed (p2, d) the
  maximizer of H1 has at most two distinct nonzero levels (KKT: ln p
  affine in p admits ≤ 2 roots; zeros recurse to lower rank), so
  enumerating two-level spectra over (rank r, multiplicity m) is
  exhaustive. The maximum coincides with the Harremoës–Topsøe family
  x = (1+√((d−1)(d·p2−1)))/d, rest equal (independent pin in tests).
- Lemma 3 (k = 3 ceiling, SOUND closed-form bound): with the size-biased
  variance σ² = p3 − p2² (σ² = 0 ⟺ flat spectrum at any rank),
  δ ≤ inf_θ [σ²/(2θ²p2²) + max_{0≤w≤η} (w·ln((d−1)p2/w) − C(θ)·w)],
  η = min(1, σ²/((1−θ)²p2²)), C(θ) = (1−θ) + (1−θ)²/(2θ²).
  Derivation: size-biased split at θp2; Taylor remainder −ln(1+u) ≤
  −u + u²/(2θ²) on the large part; Chebyshev for the small-part mass;
  concave equal-split maximum for the small-part contribution. The
  bound is 0 on every flat spectrum and continuous in (p2, p3).
- delta_max3_exact: ≤ 3-level KKT enumeration (ln p quadratic in p ⇒
  ≤ 3 interior levels). BEST-EFFORT reference for slack measurement
  ONLY — root-finding may miss branches, so it is NEVER used as the
  witness (an underestimated ceiling would be an UNSOUND witness; the
  sound object is delta_max3_bound). O27: what Lemma 3 leaves invisible
  — the gap between the bound and the true δ_max^(3); that residue is
  measured (reported), not assumed away.

Witness: U^(k) = I'2 + δ_max^(k)(SF1) + δ_max^(k)(F1F2), sound upper
bound on I(S:F2|F1) and hence on −ln F_sq(ρ, R̄(ρ_{F1F2})) (JRSWW).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from petz_witness.entropies import (
    entropy_renyi2_nats,
    entropy_vn_nats,
    i2prime_nats,
    marginals,
    purity,
    third_moment,
)
from petz_witness.linops import Matrix, Tripartition

__all__ = [
    "CONF_THETA_GRID_COUNT",
    "delta_gap",
    "delta_max2",
    "delta_max2_ht_closed_form",
    "delta_max3_bound",
    "delta_max3_exact",
    "delta_min2",
    "witness_u2",
    "witness_u3",
    "witness_u3pm",
]

#: ours: θ-grid resolution for the outer 1D infimum of Lemma 3; the
#: objective is smooth in θ and the infimum is a bound either way (any θ
#: gives a VALID bound — the grid only affects tightness, not soundness).
CONF_THETA_GRID_COUNT = 512

#: ours: moment sanity slack — p2/p3 arriving from exact simulation carry
#: O(1e-15) arithmetic noise; inputs outside [0,1] by more than this are
#: caller bugs, not noise.
_MOMENT_ATOL = 1e-12


def delta_gap(rho: Matrix) -> float:
    """Exact δ(ρ) = H1(ρ) − H2(ρ) ≥ 0 (nats)."""
    return entropy_vn_nats(rho) - entropy_renyi2_nats(rho)


def _entropy_of_spectrum(spec: NDArray[np.float64]) -> float:
    s = spec[spec > 1e-300]
    return float(-(s * np.log(s)).sum())


def delta_max2(p2: float, d: int) -> float:
    """Exact max of H1 − H2 over d-dim spectra with Tr ρ² = p2 (Lemma 2).

    Complete enumeration of two-level spectra (a×m, b×(r−m)) for all
    ranks r ≤ d. Returns 0 for p2 outside (1/d, 1] margins (flat ends).
    """
    if not (1.0 / d - _MOMENT_ATOL <= p2 <= 1.0 + _MOMENT_ATOL):
        raise ValueError(f"p2={p2} outside [1/{d}, 1]")
    p2c = min(max(p2, 1.0 / d), 1.0)
    best_h1 = -np.inf
    for r in range(1, d + 1):
        if p2c < 1.0 / r - _MOMENT_ATOL:
            continue
        if abs(p2c - 1.0 / r) <= _MOMENT_ATOL:
            best_h1 = max(best_h1, np.log(r))
            continue
        for m in range(1, r):
            a_count, b_count = float(m), float(r - m)
            # a_count·a² + (1 − a_count·a)²/b_count = p2
            c2 = a_count * b_count + a_count * a_count
            c1 = -2.0 * a_count
            c0 = 1.0 - b_count * p2c
            disc = c1 * c1 - 4.0 * c2 * c0
            if disc < 0.0:
                continue
            sq = float(np.sqrt(disc))
            for sgn in (1.0, -1.0):
                a = (-c1 + sgn * sq) / (2.0 * c2)
                b = (1.0 - a_count * a) / b_count
                if a < -_MOMENT_ATOL or b < -_MOMENT_ATOL:
                    continue
                spec = np.concatenate(
                    [np.full(m, max(a, 0.0)), np.full(r - m, max(b, 0.0))]
                )
                spec /= spec.sum()
                best_h1 = max(best_h1, _entropy_of_spectrum(spec))
    return max(best_h1 + float(np.log(p2c)), 0.0)


def delta_min2(p2: float, d: int) -> float:
    """Exact MIN of H1 − H2 over d-dim spectra with Tr ρ² = p2 (Lemma 5).

    Same complete two-level enumeration as delta_max2 (the ≤2-level KKT
    argument covers minima as well as maxima), taking the minimum.
    Registered 2026-11-13 AFTER the K2 smoke run diagnosed the one-sided
    witness's dominant cost (the discarded −δ(F1)−δ(SF1F2) cancellation)
    and BEFORE the full 10-seed map — post-smoke design extension,
    recorded per the p2-K0b precedent.
    """
    if not (1.0 / d - _MOMENT_ATOL <= p2 <= 1.0 + _MOMENT_ATOL):
        raise ValueError(f"p2={p2} outside [1/{d}, 1]")
    p2c = min(max(p2, 1.0 / d), 1.0)
    best_h1 = np.inf
    for r in range(1, d + 1):
        if p2c < 1.0 / r - _MOMENT_ATOL:
            continue
        if abs(p2c - 1.0 / r) <= _MOMENT_ATOL:
            best_h1 = min(best_h1, np.log(r))
            continue
        for m in range(1, r):
            a_count, b_count = float(m), float(r - m)
            c2 = a_count * b_count + a_count * a_count
            c1 = -2.0 * a_count
            c0 = 1.0 - b_count * p2c
            disc = c1 * c1 - 4.0 * c2 * c0
            if disc < 0.0:
                continue
            sq = float(np.sqrt(disc))
            for sgn in (1.0, -1.0):
                a = (-c1 + sgn * sq) / (2.0 * c2)
                b = (1.0 - a_count * a) / b_count
                if a < -_MOMENT_ATOL or b < -_MOMENT_ATOL:
                    continue
                spec = np.concatenate(
                    [np.full(m, max(a, 0.0)), np.full(r - m, max(b, 0.0))]
                )
                spec /= spec.sum()
                best_h1 = min(best_h1, _entropy_of_spectrum(spec))
    return max(float(best_h1) + float(np.log(p2c)), 0.0)


def delta_max2_ht_closed_form(p2: float, d: int) -> float:
    """Independent pin (M7 cross-enumeration): Harremoës–Topsøe family
    (x, (1−x)/(d−1), …), x = (1+√((d−1)(d·p2−1)))/d. Equals delta_max2 —
    tested, not assumed."""
    p2c = min(max(p2, 1.0 / d), 1.0)
    x = (1.0 + float(np.sqrt((d - 1) * (d * p2c - 1.0)))) / d
    spec = np.concatenate([[x], np.full(d - 1, (1.0 - x) / (d - 1))])
    return max(_entropy_of_spectrum(spec) + float(np.log(p2c)), 0.0)


def delta_max3_bound(p2: float, p3: float, d: int) -> float:
    """Lemma 3: sound closed-form ceiling on δ from (p2, p3, d) (nats).

    0 on every flat spectrum (σ² = 0), continuous in the moments.
    Any θ gives a valid bound; the grid infimum only tightens it.
    """
    if p2 <= 0.0:
        raise ValueError("p2 must be positive")
    sigma2 = max(p3 - p2 * p2, 0.0)
    if sigma2 == 0.0:
        return 0.0
    v = sigma2 / (p2 * p2)
    best = np.inf
    thetas = np.linspace(1.0 / (CONF_THETA_GRID_COUNT + 1), 1.0, CONF_THETA_GRID_COUNT, endpoint=False)
    a_const = (d - 1) * p2
    for theta in thetas:
        eta = min(1.0, v / ((1.0 - theta) ** 2))
        c = (1.0 - theta) + (1.0 - theta) ** 2 / (2.0 * theta * theta)
        # max over w ∈ [0, η] of g(w) = w·ln(a_const/w) − c·w (concave).
        w_star = a_const * float(np.exp(-1.0 - c))
        w_opt = min(w_star, eta)
        if w_opt <= 0.0:
            g = 0.0
        else:
            g = w_opt * (float(np.log(a_const / w_opt)) - c)
        best = min(best, v / (2.0 * theta * theta) + max(g, 0.0))
    return float(max(best, 0.0))


#: ours: residual acceptance for root candidates — scipy fsolve can
#: report ier=1 at a NON-root when the Jacobian is singular (zero step
#: satisfies the xtol criterion; reproduced 2026-11-13 on (1,2,1) at
#: p2=0.5: residual (−0.22, −0.19) returned as "converged"), so every
#: candidate is verified against the moment equations explicitly.
_ROOT_RESIDUAL_ATOL = 1e-10


def delta_max3_exact(p2: float, p3: float, d: int) -> float:
    """BEST-EFFORT exact max of δ at fixed (p2, p3, d): ≤3-level KKT
    enumeration (Lemma 4) with multi-start root finding.

    NOT sound as a witness ceiling (a missed root underestimates); used
    only to MEASURE the slack of delta_max3_bound in tests/reports.
    Every root candidate is residual-verified (see _ROOT_RESIDUAL_ATOL).
    Returns NaN when no candidate satisfies the moment constraints —
    0 would falsely read as "flat" (honest-unknown per SEMANTICS §2).
    """
    import warnings

    from scipy.optimize import fsolve

    best_h1 = -np.inf

    def consider(spec: NDArray[np.float64]) -> None:
        nonlocal best_h1
        if spec.min() < -1e-12:
            return
        spec = np.clip(spec, 0.0, None)
        if abs(spec.sum() - 1.0) > 1e-9:
            return
        if abs(float((spec**2).sum()) - p2) > 1e-9:
            return
        if abs(float((spec**3).sum()) - p3) > 1e-9:
            return
        best_h1 = max(best_h1, _entropy_of_spectrum(spec))

    # Uniform and two-level candidates (boundary of the 3-level family).
    for r in range(1, d + 1):
        consider(np.concatenate([np.full(r, 1.0 / r), np.zeros(d - r)]))
    for r in range(2, d + 1):
        for m1 in range(1, r):
            m2 = r - m1
            # two-level spectra at fixed (p1, p2); p3 checked in consider()
            a_count, b_count = float(m1), float(m2)
            c2 = a_count * b_count + a_count * a_count
            c1 = -2.0 * a_count
            c0 = 1.0 - b_count * p2
            disc = c1 * c1 - 4.0 * c2 * c0
            if disc < 0.0:
                continue
            for sgn in (1.0, -1.0):
                a = (-c1 + sgn * float(np.sqrt(disc))) / (2.0 * c2)
                b = (1.0 - a_count * a) / b_count
                consider(
                    np.concatenate([np.full(m1, a), np.full(m2, b)])
                )
    # Three-level interior candidates.
    for r in range(3, d + 1):
        for m1 in range(1, r - 1):
            for m2 in range(1, r - m1):
                m3 = r - m1 - m2
                starts = [
                    (0.6 / m1, 0.3 / m2),
                    (0.8 / m1, 0.15 / m2),
                    (0.45 / m1, 0.35 / m2),
                    (1.5 / r, 0.9 / r),
                    (2.5 / r, 0.5 / r),
                ]

                def eq3(z: NDArray[np.float64]) -> NDArray[np.float64]:
                    a, b = float(z[0]), float(z[1])
                    cc = (1.0 - m1 * a - m2 * b) / m3
                    return np.array(
                        [
                            m1 * a * a + m2 * b * b + m3 * cc * cc - p2,
                            m1 * a**3 + m2 * b**3 + m3 * cc**3 - p3,
                        ]
                    )

                for a0, b0 in starts:
                    with warnings.catch_warnings():
                        # scipy emits progress warnings on failed starts;
                        # candidates are residual-verified below, so the
                        # warning carries no information here.
                        warnings.simplefilter("ignore", RuntimeWarning)
                        sol = np.asarray(
                            fsolve(eq3, np.array([a0, b0])), dtype=np.float64
                        )
                    if float(np.abs(eq3(sol)).max()) > _ROOT_RESIDUAL_ATOL:
                        continue
                    a, b = float(sol[0]), float(sol[1])
                    cc = (1.0 - m1 * a - m2 * b) / m3
                    consider(
                        np.concatenate(
                            [np.full(m1, a), np.full(m2, b), np.full(m3, cc)]
                        )
                    )
    if not np.isfinite(best_h1):
        return float("nan")
    return max(float(best_h1) + float(np.log(p2)), 0.0)


def witness_u2(rho: Matrix, dims: Tripartition) -> float:
    """U^(2) = I'2 + δ_max^(2)(SF1) + δ_max^(2)(F1F2) from exact moments."""
    rho_sf1, rho_f1f2, _, _ = marginals(rho, dims)
    return (
        i2prime_nats(rho, dims)
        + delta_max2(purity(rho_sf1), rho_sf1.shape[0])
        + delta_max2(purity(rho_f1f2), rho_f1f2.shape[0])
    )


def witness_u3(rho: Matrix, dims: Tripartition) -> float:
    """U^(3) = I'2 + δ_max^(3)-bound(SF1) + δ_max^(3)-bound(F1F2)."""
    rho_sf1, rho_f1f2, _, _ = marginals(rho, dims)
    return (
        i2prime_nats(rho, dims)
        + delta_max3_bound(purity(rho_sf1), third_moment(rho_sf1), rho_sf1.shape[0])
        + delta_max3_bound(purity(rho_f1f2), third_moment(rho_f1f2), rho_f1f2.shape[0])
    )


def witness_u3pm(rho: Matrix, dims: Tripartition) -> float:
    """U^(3)± = U^(3) − δ_min^(2)(F1) − δ_min^(2)(SF1F2).

    Sound by the Lemma-1 identity: subtracting a certified LOWER bound
    on the discarded −δ(F1)−δ(SF1F2) pair restores part of the noise
    cancellation the one-sided witness pays for. δ_min at fixed (p2,p3)
    ≥ δ_min at fixed p2 only, so the k=2-exact minimum keeps soundness.
    """
    rho_sf1, rho_f1f2, rho_f1, rho_full = marginals(rho, dims)
    return (
        witness_u3(rho, dims)
        - delta_min2(purity(rho_f1), rho_f1.shape[0])
        - delta_min2(purity(rho_full), rho_full.shape[0])
    )
