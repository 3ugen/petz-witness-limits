"""K0 cross-check gate: two independent paths per map, CPTP on the object.

GO criterion (SKELETON K0, ratified 2026-11-13): agreement ≤ 1e-10 on
full-rank AND rank-deficient anchors; CP/TP on every constructed Choi.

M4 notes:
- power-vs-Choi agreement catches basis/embedding mistakes in either path
  (during development this gate caught an input-factor basis mismatch —
  the Choi was eigen-labeled while application was computational).
- closed-form-vs-quadrature catches a wrong spectral kernel: mutating
  k(ω) → 1 (which silently turns R̄ into R^P) fails on any non-flat
  state; flat states would NOT distinguish them — hence the non-flat
  random states here, not only anchors.
"""

from __future__ import annotations

import numpy as np
import pytest

from petz_witness import (
    QUBIT_TRIPARTITION,
    Tripartition,
    apply_f1_map_choi,
    assert_cptp,
    ghz_quantum,
    marginals,
    petz_apply_power,
    petz_choi,
    random_density,
    rho_minus,
    twirled_choi_closed_form,
    twirled_choi_quadrature,
    w_state,
)

#: SKELETON K0 GO threshold (ratified 2026-11-13), dimensionless.
CONF_XCHECK_ATOL = 1e-10
#: ours: quadrature path tolerance — limited by Gauss–Legendre truncation,
#: see CONF_QUAD_* in petz.py; stricter than any downstream use.
CONF_QUAD_ATOL = 1e-8

DIMS = QUBIT_TRIPARTITION


def _random_tripartite(seed: int, rank: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return random_density(DIMS.d_total, rank, rng)


CASES: list[tuple[str, np.ndarray]] = [
    ("rho_minus(rank-deficient)", rho_minus()),
    ("ghz(rank-1)", ghz_quantum()),
    ("w(rank-1)", w_state()),
    ("random-full", _random_tripartite(11, 8)),
    ("random-rank3", _random_tripartite(12, 3)),
    ("random-full-2", _random_tripartite(13, 8)),
]


@pytest.mark.parametrize(("name", "rho"), CASES, ids=[c[0] for c in CASES])
def test_petz_power_vs_choi(name: str, rho: np.ndarray) -> None:
    rho_sf1, rho_f1f2, rho_f1, _ = marginals(rho, DIMS)
    out_power = petz_apply_power(rho_sf1, rho_f1, rho_f1f2, DIMS)
    choi = petz_choi(rho_sf1, rho_f1, DIMS)
    assert_cptp(choi, DIMS, atol=CONF_XCHECK_ATOL)
    out_choi = apply_f1_map_choi(choi, rho_f1f2, DIMS)
    assert float(np.abs(out_power - out_choi).max()) < CONF_XCHECK_ATOL


@pytest.mark.parametrize(("name", "rho"), CASES, ids=[c[0] for c in CASES])
def test_twirled_closed_form_vs_quadrature(name: str, rho: np.ndarray) -> None:
    rho_sf1, _, rho_f1, _ = marginals(rho, DIMS)
    choi_cf = twirled_choi_closed_form(rho_sf1, rho_f1, DIMS)
    choi_q = twirled_choi_quadrature(rho_sf1, rho_f1, DIMS)
    assert_cptp(choi_cf, DIMS, atol=1e-9)
    assert float(np.abs(choi_cf - choi_q).max()) < CONF_QUAD_ATOL


def test_twirled_equals_petz_on_flat_manifold() -> None:
    """On flat ρ_{SF1}, ρ_{F1} all rotations are scalar phases: R̄ = R^P."""
    rho = ghz_quantum()
    rho_sf1, _, rho_f1, _ = marginals(rho, DIMS)
    choi_tw = twirled_choi_closed_form(rho_sf1, rho_f1, DIMS)
    choi_p = petz_choi(rho_sf1, rho_f1, DIMS)
    assert float(np.abs(choi_tw - choi_p).max()) < CONF_XCHECK_ATOL


def test_twirled_differs_from_petz_off_flat_manifold() -> None:
    """Positive control for the kernel (M4): on a non-flat random state the
    two maps must DIFFER — a k(ω) ≡ 1 mutation would pass the flat test
    and fail here."""
    rho = _random_tripartite(21, 8)
    rho_sf1, _, rho_f1, _ = marginals(rho, DIMS)
    choi_tw = twirled_choi_closed_form(rho_sf1, rho_f1, DIMS)
    choi_p = petz_choi(rho_sf1, rho_f1, DIMS)
    assert float(np.abs(choi_tw - choi_p).max()) > 1e-4


def test_asymmetric_dimensions() -> None:
    """Non-square tripartition (d_s, d_1, d_2) = (2, 4, 2): embedding and
    reshape orders are exercised beyond the all-qubit symmetric case —
    a transposed kron order passes 2×2×2 by accident and fails here."""
    dims = Tripartition(2, 4, 2)
    rng = np.random.default_rng(31)
    rho = random_density(dims.d_total, dims.d_total, rng)
    from petz_witness.entropies import marginals as marg

    rho_sf1, rho_f1f2, rho_f1, _ = marg(rho, dims)
    out_power = petz_apply_power(rho_sf1, rho_f1, rho_f1f2, dims)
    choi = petz_choi(rho_sf1, rho_f1, dims)
    assert_cptp(choi, dims, atol=CONF_XCHECK_ATOL)
    out_choi = apply_f1_map_choi(choi, rho_f1f2, dims)
    assert float(np.abs(out_power - out_choi).max()) < CONF_XCHECK_ATOL
