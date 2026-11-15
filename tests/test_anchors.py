"""Anchor tests: every expected value comes from a source OTHER than the
code under test (R-PINNED-BOUND): hand derivations in
notes/2026-11-13-carbon-v5-p3-skeleton-analysis.md and closed forms.

M4 notes (wrong implementation each test would catch):
- rho⁻/rho⁺ I'2 signs → an implementation that clamps I'2 at 0, or swaps
  marginals, fails here.
- Bell anchor −ln F_sq = CMI → a ROOT-fidelity implementation (reports
  −2 ln F = 2·CMI) fails here; so does a log2-based one (factor ln 2).
- Petz-exactness on classical MCs → a map built with the wrong embedding
  order (F1 vs F2) fails here.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petz_witness import (
    QUBIT_TRIPARTITION,
    bell_pair_trivial_f1,
    cmi_nats,
    fidelity_sq,
    ghz_classical,
    ghz_quantum,
    i2prime_nats,
    marginals,
    neg_log_fidelity_sq,
    petz_apply_power,
    product_state,
    rho_minus,
    rho_plus,
    twirled_choi_closed_form,
    apply_f1_map_choi,
    w_state,
)

DIMS = QUBIT_TRIPARTITION
ATOL = 1e-10


def _petz_recovered(rho: np.ndarray) -> np.ndarray:
    rho_sf1, rho_f1f2, rho_f1, _ = marginals(rho, DIMS)
    return petz_apply_power(rho_sf1, rho_f1, rho_f1f2, DIMS)


def _twirled_recovered(rho: np.ndarray) -> np.ndarray:
    rho_sf1, rho_f1f2, rho_f1, _ = marginals(rho, DIMS)
    choi = twirled_choi_closed_form(rho_sf1, rho_f1, DIMS)
    return apply_f1_map_choi(choi, rho_f1f2, DIMS)


@pytest.mark.parametrize(
    ("state", "expected_i2p"),
    [
        (rho_minus, math.log(8.0 / 9.0)),
        (rho_plus, math.log(10.0 / 9.0)),
    ],
)
def test_classical_markov_chains_i2prime_both_signs(state, expected_i2p) -> None:
    """Carbon v5 §4 counterexamples: CMI = 0, Petz-exact, I'2 ≠ 0."""
    rho = state()
    assert abs(cmi_nats(rho, DIMS)) < ATOL
    assert i2prime_nats(rho, DIMS) == pytest.approx(expected_i2p, abs=1e-12)
    rec = _petz_recovered(rho)
    assert fidelity_sq(rho, rec) == pytest.approx(1.0, abs=1e-9)
    rec_tw = _twirled_recovered(rho)
    assert fidelity_sq(rho, rec_tw) == pytest.approx(1.0, abs=1e-9)


def test_ghz_quantum_flat_manifold() -> None:
    """GHZ: CMI = ln 2 = I'2 (flat marginals ⇒ surrogate exact)."""
    rho = ghz_quantum()
    assert cmi_nats(rho, DIMS) == pytest.approx(math.log(2.0), abs=1e-12)
    assert i2prime_nats(rho, DIMS) == pytest.approx(math.log(2.0), abs=1e-12)


def test_ghz_classical_is_markov() -> None:
    rho = ghz_classical()
    assert abs(cmi_nats(rho, DIMS)) < ATOL
    assert abs(i2prime_nats(rho, DIMS)) < ATOL
    assert fidelity_sq(rho, _petz_recovered(rho)) == pytest.approx(1.0, abs=1e-9)


def test_w_state_cmi_closed_form() -> None:
    """W₃: CMI = ln 3 − (2/3)ln 2 from marginal spectra {2/3, 1/3}."""
    expected = math.log(3.0) - (2.0 / 3.0) * math.log(2.0)
    assert cmi_nats(w_state(), DIMS) == pytest.approx(expected, abs=1e-12)


def test_product_state_trivial() -> None:
    rho = product_state()
    assert abs(cmi_nats(rho, DIMS)) < ATOL
    assert abs(i2prime_nats(rho, DIMS)) < ATOL
    assert fidelity_sq(rho, _petz_recovered(rho)) == pytest.approx(1.0, abs=1e-9)


def test_bell_anchor_pins_fidelity_convention() -> None:
    """−ln F_sq(ρ, R̄(ρ_{F1F2})) = CMI = 2 ln 2 exactly on the Bell anchor.

    The recovery of a trivial-F1 state is ρ_{SF1} ⊗-structured guesswork:
    R̄(ρ_{F1F2}) = ρ_{SF1} ⊗ ρ_{F2}-recovery = maximally mixed on S F2
    side, giving F_sq = 1/4. A root-fidelity implementation reports
    2·CMI here and fails (M4 wrong-implementation).
    """
    rho = bell_pair_trivial_f1()
    cmi = cmi_nats(rho, DIMS)
    assert cmi == pytest.approx(2.0 * math.log(2.0), abs=1e-12)
    rec = _twirled_recovered(rho)
    assert neg_log_fidelity_sq(rho, rec) == pytest.approx(cmi, abs=1e-9)
