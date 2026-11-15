"""K2 ensemble sanity gates: properties VERIFIED against the exact
oracles, not assumed from the literature (T4).

M4 notes: the HJPW exactness test fails if the block structure is built
on the wrong factor (F1R instead of F1L); the cluster-chain QMC test
fails if the CZ chain or the contiguous partition is mis-indexed; the
noise test pins the named channel (a global-depolarizing implementation
gives a different single-qubit marginal and fails).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petz_witness import cmi_nats, fidelity_sq, marginals, petz_apply_power
from petz_witness.ensembles import (
    all_ensembles,
    apply_per_qubit_depolarizing,
    ensemble_a_stabilizer_perturbed,
    ensemble_f_hjpw,
)
from petz_witness.linops import Tripartition, partial_trace


def test_seeded_reproducibility() -> None:
    e1 = all_ensembles(7)
    e2 = all_ensembles(7)
    for k in e1:
        assert len(e1[k]) == len(e2[k])
        for s1, s2 in zip(e1[k], e2[k]):
            assert s1.label == s2.label
            assert np.array_equal(s1.rho, s2.rho)


def test_states_are_density_matrices() -> None:
    for states in all_ensembles(3).values():
        for st in states:
            rho = st.rho
            assert rho.shape == (st.dims.d_total, st.dims.d_total)
            assert np.trace(rho).real == pytest.approx(1.0, abs=1e-9)
            w = np.linalg.eigvalsh(rho)
            assert float(w.min()) > -1e-10


def test_cluster_chain_is_qmc() -> None:
    """Verified, not assumed: 1D cluster states have I(S:F2|F1) = 0 for
    contiguous tripartitions and are Petz-recoverable exactly."""
    checked = 0
    for st in ensemble_a_stabilizer_perturbed(0):
        if "cluster" not in st.label or ":eps0.0:" not in st.label:
            continue
        checked += 1
        assert abs(cmi_nats(st.rho, st.dims)) < 1e-9, st.label
        rho_sf1, rho_f1f2, rho_f1, _ = marginals(st.rho, st.dims)
        rec = petz_apply_power(rho_sf1, rho_f1, rho_f1f2, st.dims)
        # tolerance 2e-6: fidelity of rank-deficient states carries
        # √(eigh noise) ≈ 1e-7 per near-zero eigenvalue at d = 32
        # (sqrt amplification at the spectrum edge, not a map error)
        assert fidelity_sq(st.rho, rec) == pytest.approx(1.0, abs=2e-6), st.label
    assert checked == 2  # cluster-4 and cluster-5 unperturbed bases


def test_entangling_perturbation_moves_cmi() -> None:
    """The ε > 0 states must actually leave the Markov manifold —
    product rotations would keep CMI = 0 (the M4 wrong-implementation
    observed on the first draft of this ensemble)."""
    moved = 0
    for st in ensemble_a_stabilizer_perturbed(0):
        if ":eps0.2:" in st.label:
            if abs(cmi_nats(st.rho, st.dims)) > 1e-3:
                moved += 1
    assert moved > 0


def test_hjpw_states_are_exact_qmc_with_nontrivial_mi() -> None:
    for st in ensemble_f_hjpw(11):
        assert abs(cmi_nats(st.rho, st.dims)) < 1e-9, st.label
        rho_sf1, rho_f1f2, rho_f1, _ = marginals(st.rho, st.dims)
        rec = petz_apply_power(rho_sf1, rho_f1, rho_f1f2, st.dims)
        assert fidelity_sq(st.rho, rec) == pytest.approx(1.0, abs=1e-7), st.label
        d = (st.dims.d_s, st.dims.d_1, st.dims.d_2)
        from petz_witness.entropies import entropy_vn_nats

        mi = (
            entropy_vn_nats(partial_trace(st.rho, d, (0,)))
            + entropy_vn_nats(partial_trace(st.rho, d, (1,)))
            - entropy_vn_nats(partial_trace(st.rho, d, (0, 1)))
        )
        assert mi >= 0.5 * math.log(2.0) - 1e-9, st.label


def test_noise_channel_pins_single_qubit_marginal() -> None:
    """Per-qubit depolarizing: |0⟩⟨0| marginal → diag(1−p/2, p/2)."""
    dims = Tripartition(2, 2, 2)
    rho = np.zeros((8, 8), dtype=np.complex128)
    rho[0, 0] = 1.0
    p = 0.12
    noisy = apply_per_qubit_depolarizing(rho, dims, p)
    rho_s = partial_trace(noisy, (2, 2, 2), (0,))
    assert rho_s[0, 0].real == pytest.approx(1.0 - p / 2.0, abs=1e-12)
    assert rho_s[1, 1].real == pytest.approx(p / 2.0, abs=1e-12)
