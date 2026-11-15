"""K0 theorem-sanity gate: I(S:F2|F1) ≥ −ln F_sq(ρ, R̄(ρ_{F1F2})) − 1e-12
on 500 seeded random states (SKELETON K0 GO, ratified 2026-11-13).

Source of the inequality (O22): JRSWW 1509.07127 Theorem 2.1 with
N = Tr_S, in the F_sq convention (−ln F_sq, NOT −2 ln F_sq — the Bell
anchor in test_anchors.py pins the factor).

M7 note: 500 random 3-qubit states do not PROVE the theorem (it is
proved in the source); this gate detects implementation bugs.
M4 calibration (mutations actually run, 2026-11-13):
- root-fidelity convention slip (margin computed as CMI − 2·(−ln F_sq)):
  RED on 100/100 random states — this gate catches it;
- the half-log kernel bug (ω/2)/sinh(ω/2): SURVIVOR here (worst margin
  +0.079 — the buggy map is still a channel near Petz and does not
  violate the bound empirically). It is NOT undistinguished overall:
  test_twirled_closed_form_vs_quadrature went RED on exactly this bug
  in the real development run; that gate, not this one, owns kernel
  correctness. Division of labor recorded per M4.
"""

from __future__ import annotations

import numpy as np

from petz_witness import (
    QUBIT_TRIPARTITION,
    apply_f1_map_choi,
    cmi_nats,
    marginals,
    neg_log_fidelity_sq,
    random_density,
    twirled_choi_closed_form,
)

#: SKELETON K0 (ratified 2026-11-13): sample size of the sanity gate.
CONF_SANITY_STATES_COUNT = 500
#: SKELETON K0 (ratified 2026-11-13): allowed numerical slack, nats.
CONF_SANITY_SLACK = 1e-12
#: ours: master seed of the sanity ensemble — fixed so the gate is
#: reproducible; independence from the implementation comes from the
#: theorem, not the seed.
CONF_SANITY_SEED = 20261113

DIMS = QUBIT_TRIPARTITION


def test_jrsww_bound_on_random_states() -> None:
    rng = np.random.default_rng(CONF_SANITY_SEED)
    worst_margin = float("inf")
    for _ in range(CONF_SANITY_STATES_COUNT):
        rank = int(rng.integers(1, DIMS.d_total + 1))
        rho = random_density(DIMS.d_total, rank, rng)
        rho_sf1, rho_f1f2, rho_f1, _ = marginals(rho, DIMS)
        choi = twirled_choi_closed_form(rho_sf1, rho_f1, DIMS)
        rec = apply_f1_map_choi(choi, rho_f1f2, DIMS)
        margin = cmi_nats(rho, DIMS) - neg_log_fidelity_sq(rho, rec)
        worst_margin = min(worst_margin, margin)
    assert worst_margin >= -CONF_SANITY_SLACK, (
        f"JRSWW bound violated: worst margin {worst_margin:.3e} nats"
    )
