"""K3 estimator gates: unbiasedness against exact moments (T4 — the
exact oracles are the independent reference), factor-table pins, and
shadow-reconstruction sanity.

M4 notes: the pair-table pin fails if the shadow inverse channel is
mis-scaled (e.g., 2P−I instead of 3P−I gives Tr(σσ') = 1 on equal
outcomes, not 5); the unbiasedness gate fails on an off-by-one in the
U-statistic denominators (diagonal pairs included → p̂2 biased high by
~Tr(σ²)/N); the probability-table test fails if the qubit/digit order
in basis decoding is transposed.
"""

from __future__ import annotations

import numpy as np
import pytest

from petz_witness import QUBIT_TRIPARTITION, ghz_classical, marginals
from petz_witness.entropies import purity, third_moment
from petz_witness.linops import Tripartition, random_density
from petz_witness.sampling import (
    _PAIR_TABLE,
    _probability_table,
    estimate_p2_blocks,
    estimate_p3_blocks,
    sample_dataset,
    shadow_reconstruction,
    witness_ucb_from_dataset,
)

DIMS = QUBIT_TRIPARTITION


def test_pair_table_pinned_values() -> None:
    """Same basis+bit → 5; same basis, flipped bit → −4; cross → 1/2.
    (Pinned by hand from Tr(σσ') = 9|⟨ψ|ψ'⟩|² − 4.)"""
    assert _PAIR_TABLE[0, 0] == pytest.approx(5.0, abs=1e-12)
    assert _PAIR_TABLE[0, 1] == pytest.approx(-4.0, abs=1e-12)
    assert _PAIR_TABLE[0, 2] == pytest.approx(0.5, abs=1e-12)
    assert _PAIR_TABLE[4, 1] == pytest.approx(0.5, abs=1e-12)


def test_probability_table_matches_exact_marginals() -> None:
    """All-Z basis row of the table = diagonal of ρ (classical GHZ)."""
    rho = ghz_classical()
    table = _probability_table(rho, 3)
    z_index = 3**3 - 1  # digits (2,2,2) = Z on every qubit
    np.testing.assert_allclose(
        table[z_index], np.real(np.diag(rho)), atol=1e-12
    )


def test_moment_estimators_unbiased() -> None:
    """p̂2/p̂3 of every marginal converge to exact moments (3-qubit
    random mixed state, N = 60k): within 4 block-SEs."""
    rng = np.random.default_rng(90)
    rho = random_density(8, 5, rng)
    ds = sample_dataset(rho, DIMS, 60_000, rng)
    rho_sf1, rho_f1f2, rho_f1, rho_full = marginals(rho, DIMS)
    cases = [
        ((0, 1), rho_sf1),
        ((1, 2), rho_f1f2),
        ((1,), rho_f1),
        ((0, 1, 2), rho_full),
    ]
    for qubits, marg in cases:
        blocks = estimate_p2_blocks(ds, qubits)
        est, se = blocks.mean(), blocks.std(ddof=1) / np.sqrt(len(blocks))
        assert est == pytest.approx(purity(marg), abs=max(4 * se, 1e-3)), qubits
    for qubits, marg in cases[:2]:
        blocks = estimate_p3_blocks(ds, qubits, rng)
        est, se = blocks.mean(), blocks.std(ddof=1) / np.sqrt(len(blocks))
        assert est == pytest.approx(third_moment(marg), abs=max(4 * se, 2e-3)), qubits


def test_witness_ucb_covers_exact_u3pm() -> None:
    """One-sided 95% UCB covers the exact U^(3)± (checked on a certified-
    region state where the witness value matters)."""
    from petz_witness.certificate import witness_u3pm
    from petz_witness.ensembles import apply_per_qubit_depolarizing

    rng = np.random.default_rng(91)
    rho = apply_per_qubit_depolarizing(ghz_classical(), DIMS, 1e-3)
    exact = witness_u3pm(rho, DIMS)
    covered = 0
    for rep in range(10):
        ds = sample_dataset(rho, DIMS, 8_000, rng)
        _, ucb = witness_ucb_from_dataset(ds, DIMS, rng)
        if ucb >= exact - 1e-9:
            covered += 1
    assert covered >= 8  # nominal 95% one-sided; 10 reps, allow 2 misses


def test_shadow_reconstruction_converges() -> None:
    rng = np.random.default_rng(92)
    rho = random_density(8, 3, rng)
    ds = sample_dataset(rho, DIMS, 40_000, rng)
    rho_hat = shadow_reconstruction(ds, DIMS)
    assert float(np.abs(rho_hat - rho).max()) < 0.05
    w = np.linalg.eigvalsh(rho_hat)
    assert float(w.min()) >= -1e-12
    assert float(np.trace(rho_hat).real) == pytest.approx(1.0, abs=1e-9)


def test_newton_girard_exact_vs_brute_force() -> None:
    """v4 NG complete U-statistics equal brute-force sums (probe B of
    the v6 analysis, kept as a gate; M4: a wrong diagonal-correction
    constant (5^k/7^k) shifts both moments and fails at 1e-12)."""
    from petz_witness.sampling import _PAIR_TABLE, _TRIPLE_TABLE, _marginal_ng

    rng = np.random.default_rng(3)
    rho = random_density(8, 4, rng)
    ds = sample_dataset(rho, DIMS, 50, rng)
    codes, qubits = ds.codes, (0, 1)
    n = codes.shape[0]
    p2_bf = 0.0
    p3_bf = 0.0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            pp = 1.0
            for q in qubits:
                pp *= _PAIR_TABLE[codes[i, q], codes[j, q]]
            p2_bf += pp
            for k in range(n):
                if k == i or k == j:
                    continue
                pc = 1.0 + 0.0j
                for q in qubits:
                    pc *= _TRIPLE_TABLE[codes[i, q], codes[j, q], codes[k, q]]
                p3_bf += pc.real
    p2_bf /= n * (n - 1)
    p3_bf /= n * (n - 1) * (n - 2)
    m = _marginal_ng(codes, qubits)
    assert m.p2 == pytest.approx(p2_bf, abs=1e-12)
    assert m.p3 == pytest.approx(p3_bf, abs=1e-12)


def test_v4_ucb_covers_and_beats_v3() -> None:
    """v4 UCB covers exact U^(3)± and is far tighter than v3 on a
    certified-region state (the Monte-Carlo floor is gone)."""
    from petz_witness.certificate import witness_u3pm
    from petz_witness.ensembles import _classical_ghz
    from petz_witness.sampling import witness_ucb_v4

    rng = np.random.default_rng(44)
    rho = _classical_ghz(3)
    exact = witness_u3pm(rho, DIMS)
    covered = 0
    for _ in range(6):
        ds = sample_dataset(rho, DIMS, 20_000, rng)
        _, ucb = witness_ucb_v4(ds, DIMS)
        if ucb >= exact - 1e-9:
            covered += 1
        assert ucb < 1.0  # v3's floor-driven UCBs were ~2 nats
    assert covered >= 5
