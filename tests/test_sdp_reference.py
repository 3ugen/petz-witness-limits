"""Reference oracle: CPTP-SDP fidelity of recovery ≥ twirled-Petz fidelity.

The SDP optimizes over ALL CPTP recoveries, so on every state
F_FoR ≥ F_sq(ρ, R̄(ρ_{F1F2})) up to solver accuracy; on exact Markov
chains both equal 1. A violation beyond solver tolerance means a bug in
one of the two (M4: this pair is a mutual gate — e.g. a non-TP Choi in
the SDP constraint set inflates F_FoR above 1 and fails the cap check).
"""

from __future__ import annotations

import numpy as np
import pytest

from petz_witness import (
    QUBIT_TRIPARTITION,
    apply_f1_map_choi,
    fidelity_sq,
    ghz_quantum,
    marginals,
    random_density,
    rho_minus,
    twirled_choi_closed_form,
)
from petz_witness.sdp import fidelity_of_recovery_sq

#: ours: SCS first-order solver accuracy on these problem sizes (observed
#: during calibration; conservative). Dimensionless fidelity units.
CONF_SDP_ATOL = 5e-4

DIMS = QUBIT_TRIPARTITION


def _twirled_fid(rho: np.ndarray) -> float:
    rho_sf1, rho_f1f2, rho_f1, _ = marginals(rho, DIMS)
    choi = twirled_choi_closed_form(rho_sf1, rho_f1, DIMS)
    return fidelity_sq(rho, apply_f1_map_choi(choi, rho_f1f2, DIMS))


def test_sdp_equals_one_on_markov_chain() -> None:
    f = fidelity_of_recovery_sq(rho_minus(), DIMS)
    assert f == pytest.approx(1.0, abs=CONF_SDP_ATOL)


def test_sdp_dominates_twirled_petz() -> None:
    states = [ghz_quantum(), random_density(8, 8, np.random.default_rng(41))]
    for rho in states:
        f_sdp = fidelity_of_recovery_sq(rho, DIMS)
        f_tw = _twirled_fid(rho)
        assert f_sdp >= f_tw - CONF_SDP_ATOL
        assert f_sdp <= 1.0 + CONF_SDP_ATOL
