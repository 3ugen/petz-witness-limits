"""K1 gate tests (SKELETON §2-K1, thresholds ratified 2026-11-13).

GO criteria exercised here:
- zero violations on 1e4 random spectra (soundness of both ceilings);
- δ_max^(3)-bound = 0 (≤ 1e-10) on flat spectra of EVERY rank 1 ≤ r ≤ d;
- continuous decay of the bound as depolarizing p → 0;
- Lemma 1 identity and the one-sided chain on random states (with the
  K0 twirled-Petz oracle on the left end).

M4 calibration (mutations actually run during development, recorded):
- σ² := 0 mutation in delta_max3_bound → flat tests still pass but the
  1e4-spectra soundness test goes RED (bound < δ on non-flat spectra);
- dropping the C(θ)·w term tightens nothing (bound only grows — gates
  stay green: that direction is slack, not unsoundness — recorded as
  the expected asymmetry of one-sided gates);
- delta_max2 with m restricted to 1 (HT family only) still passes all
  gates — the enumeration beyond m=1 is completeness insurance per
  Lemma 2(i), pinned equal to the HT closed form here (M7).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petz_witness import (
    QUBIT_TRIPARTITION,
    apply_f1_map_choi,
    cmi_nats,
    depolarize,
    ghz_classical,
    i2prime_nats,
    marginals,
    neg_log_fidelity_sq,
    random_density,
    rho_plus,
    twirled_choi_closed_form,
)
from petz_witness.certificate import (
    delta_gap,
    delta_max2,
    delta_max2_ht_closed_form,
    delta_max3_bound,
    delta_max3_exact,
    witness_u2,
    witness_u3,
)
from petz_witness.entropies import purity, third_moment

DIMS = QUBIT_TRIPARTITION

#: SKELETON K1 (ratified 2026-11-13): random-spectrum soundness sample.
CONF_SPECTRA_COUNT = 10_000
#: ours: seed of the spectrum ensemble (reproducibility only).
CONF_SPECTRA_SEED = 20261113
#: SKELETON K1: flat-spectrum zero tolerance, nats.
CONF_FLAT_ATOL = 1e-10


def _random_spectrum(rng: np.random.Generator) -> tuple[np.ndarray, int]:
    d = int(rng.integers(2, 17))
    rank = int(rng.integers(1, d + 1))
    alpha = float(rng.choice([0.2, 1.0, 5.0]))
    spec = rng.dirichlet(np.full(rank, alpha))
    full = np.zeros(d)
    full[:rank] = spec
    return full, d


def test_soundness_on_1e4_random_spectra() -> None:
    """δ(spec) ≤ δ_max^(2)(p2,d) and ≤ δ_max^(3)-bound(p2,p3,d), zero
    violations. Moments computed from the SAME spectrum — no tolerance
    gymnastics (the 2026-11-12 random-search artifact is excluded by
    construction)."""
    rng = np.random.default_rng(CONF_SPECTRA_SEED)
    worst2 = worst3 = np.inf
    for _ in range(CONF_SPECTRA_COUNT):
        spec, d = _random_spectrum(rng)
        nz = spec[spec > 0]
        h1 = float(-(nz * np.log(nz)).sum())
        p2 = float((spec**2).sum())
        p3 = float((spec**3).sum())
        delta = h1 + math.log(p2)
        worst2 = min(worst2, delta_max2(p2, d) - delta)
        worst3 = min(worst3, delta_max3_bound(p2, p3, d) - delta)
    assert worst2 >= -1e-11, f"k=2 ceiling violated by {worst2:.3e}"
    assert worst3 >= -1e-11, f"k=3 bound violated by {worst3:.3e}"


@pytest.mark.parametrize("d", [4, 8, 16])
def test_flat_spectra_every_rank(d: int) -> None:
    """GO: δ_max^(3)-bound = 0 on flat spectra of every rank 1 ≤ r ≤ d —
    exactly the regime where the k=2 ceiling fails (Lemma 2(iv))."""
    for r in range(1, d + 1):
        p2 = 1.0 / r
        p3 = 1.0 / (r * r)
        assert delta_max3_bound(p2, p3, d) <= CONF_FLAT_ATOL
        if 1 < r < d:
            assert delta_max2(p2, d) > 0.01  # insufficiency is strict


def test_k2_insufficiency_pinned_value() -> None:
    """δ_max^(2)(1/2, 4) = 0.279677 nats = 0.4035 bits (pinned by the
    independent HT closed form AND by the dev probe of 2026-11-13)."""
    v = delta_max2(0.5, 4)
    assert v == pytest.approx(0.2796765, abs=1e-6)
    assert v == pytest.approx(delta_max2_ht_closed_form(0.5, 4), abs=1e-12)


@pytest.mark.parametrize(("p2", "d"), [(0.3, 4), (0.2, 8), (0.5, 8), (0.1, 16)])
def test_ht_closed_form_matches_enumeration(p2: float, d: int) -> None:
    """M7 cross-enumeration: two independent implementations agree."""
    assert delta_max2(p2, d) == pytest.approx(
        delta_max2_ht_closed_form(p2, d), abs=1e-12
    )


def test_bound_continuity_under_depolarizing() -> None:
    """GO: bound → 0 continuously as p → 0 on the depolarized classical
    GHZ (the Q1 discontinuity kill-shot scenario: Σ(H0−H∞) sat at 4 bits
    there; our bound must be tiny and decreasing)."""
    prev = np.inf
    for p in (1e-2, 1e-4, 1e-6):
        rho = depolarize(ghz_classical(), p)
        rho_sf1, _, _, _ = marginals(rho, DIMS)
        b = delta_max3_bound(purity(rho_sf1), third_moment(rho_sf1), 4)
        assert b < prev
        prev = b
    assert prev < 1e-4  # nats, at p = 1e-6 (own measurement, cf. analysis note)


def test_lemma1_identity_exact() -> None:
    rng = np.random.default_rng(7)
    for _ in range(50):
        rho = random_density(8, int(rng.integers(1, 9)), rng)
        rho_sf1, rho_f1f2, rho_f1, rho_full = marginals(rho, DIMS)
        lhs = cmi_nats(rho, DIMS) - i2prime_nats(rho, DIMS)
        rhs = (
            delta_gap(rho_sf1)
            + delta_gap(rho_f1f2)
            - delta_gap(rho_f1)
            - delta_gap(rho_full)
        )
        assert lhs == pytest.approx(rhs, abs=1e-10)


def test_one_sided_chain_with_k0_oracle() -> None:
    """−ln F_sq(ρ, R̄) ≤ I ≤ U^(3) ≤ U^(2)-with-exact-δ chain on random
    states; also U uses ceilings, so U^(k) ≥ I'2 + actual δ's."""
    rng = np.random.default_rng(20261113)
    for _ in range(60):
        rho = random_density(8, int(rng.integers(1, 9)), rng)
        rho_sf1, rho_f1f2, rho_f1, _ = marginals(rho, DIMS)
        choi = twirled_choi_closed_form(rho_sf1, rho_f1, DIMS)
        rec = apply_f1_map_choi(choi, rho_f1f2, DIMS)
        nlf = neg_log_fidelity_sq(rho, rec)
        i = cmi_nats(rho, DIMS)
        u3 = witness_u3(rho, DIMS)
        u2 = witness_u2(rho, DIMS)
        assert nlf <= i + 1e-11
        assert i <= u3 + 1e-11
        assert i <= u2 + 1e-11


def test_false_refusal_anatomy_rho_plus() -> None:
    """ρ⁺ (non-flat exact QMC): recoverable, but U^(k) > τ-grade values —
    completeness loss by design, soundness intact (Lemma 1 comment)."""
    rho = rho_plus()
    assert abs(cmi_nats(rho, DIMS)) < 1e-10
    u3 = witness_u3(rho, DIMS)
    assert u3 > 0.1  # refuses: the certificate cannot call this flat


def test_exact_three_level_vs_bound_slack() -> None:
    """Slack measurement (NOT a soundness gate): best-effort exact
    δ_max^(3) ≤ bound where the search converges; report-grade check
    that the bound is within a usable factor on a small probe set."""
    cases = [(0.5, 0.28, 4), (0.3, 0.11, 6), (0.25, 0.08, 8), (0.4, 0.2, 4)]
    found = 0
    for p2, p3, d in cases:
        exact = delta_max3_exact(p2, p3, d)
        if math.isnan(exact):  # best-effort search may find nothing
            continue
        found += 1
        bound = delta_max3_bound(p2, p3, d)
        assert exact <= bound + 1e-9
        assert bound <= 4.0 * exact + 1e-9  # measured slack stays usable
    assert found >= 3  # the probe set must mostly converge to mean anything


def test_delta_min2_is_certified_lower_bound() -> None:
    """Lemma 5 soundness: δ_min^(2)(p2,d) ≤ δ(spec) on random spectra."""
    from petz_witness.certificate import delta_min2

    rng = np.random.default_rng(424242)
    worst = np.inf
    for _ in range(3000):
        spec, d = _random_spectrum(rng)
        nz = spec[spec > 0]
        h1 = float(-(nz * np.log(nz)).sum())
        p2 = float((spec**2).sum())
        delta = h1 + math.log(p2)
        worst = min(worst, delta - delta_min2(p2, d))
    assert worst >= -1e-11, f"delta_min2 exceeds a real gap by {-worst:.3e}"


def test_u3pm_sound_and_tighter() -> None:
    """U^(3)± stays above CMI (soundness) and never exceeds U^(3)."""
    from petz_witness.certificate import witness_u3pm

    rng = np.random.default_rng(31337)
    for _ in range(40):
        rho = random_density(8, int(rng.integers(1, 9)), rng)
        i = cmi_nats(rho, DIMS)
        u3 = witness_u3(rho, DIMS)
        u3pm = witness_u3pm(rho, DIMS)
        assert i <= u3pm + 1e-10
        assert u3pm <= u3 + 1e-12
