"""Entropic quantities for the tripartition S|F1|F2 — all in NATS.

Sources (O22): CMI decomposition I(S:F2|F1) = S(SF1)+S(F1F2)−S(F1)−S(SF1F2)
is standard (e.g. Sutter, "Approximate Quantum Markov Chains", 1907.09421
eq. 2.5); the surrogate I'2 is the SAME linear combination with Rényi-2
entropies — explicitly NOT a Rényi generalization of CMI (BSW 1403.6102
eq. 1.3 and the discussion rejecting such formulas; LMW 1212.0248: no
universal sign or monotonicity). Consumers: K1/K2 gates of SKELETON.md,
which rely on exactly this sign structure via the identity
I − I'2 = δ(SF1)+δ(F1F2)−δ(F1)−δ(SF1F2), δ(M) = H1(M)−H2(M) ≥ 0.
"""

from __future__ import annotations

import numpy as np

from petz_witness.linops import Matrix, Tripartition, func_hermitian, partial_trace

__all__ = [
    "entropy_vn_nats",
    "entropy_renyi2_nats",
    "purity",
    "third_moment",
    "marginals",
    "cmi_nats",
    "i2prime_nats",
    "sqrtm_psd",
]


def entropy_vn_nats(rho: Matrix) -> float:
    """Von Neumann entropy −Tr ρ ln ρ in nats (0·ln 0 := 0)."""
    w = np.linalg.eigvalsh(rho)
    w = np.clip(w, 0.0, None)
    nz = w[w > 0.0]
    return float(-(nz * np.log(nz)).sum())


def purity(rho: Matrix) -> float:
    """p2 = Tr ρ² (real part; ρ Hermitian)."""
    return float(np.real(np.trace(rho @ rho)))


def third_moment(rho: Matrix) -> float:
    """p3 = Tr ρ³ — the k=3 certificate input (σ² = p3 − p2², K1)."""
    return float(np.real(np.trace(rho @ rho @ rho)))


def entropy_renyi2_nats(rho: Matrix) -> float:
    """H2 = −ln Tr ρ² in nats."""
    return float(-np.log(purity(rho)))


def marginals(rho: Matrix, dims: Tripartition) -> tuple[Matrix, Matrix, Matrix, Matrix]:
    """(ρ_SF1, ρ_F1F2, ρ_F1, ρ_SF1F2) — the four marginals of the witness."""
    d = (dims.d_s, dims.d_1, dims.d_2)
    rho_sf1 = partial_trace(rho, d, (0, 1))
    rho_f1f2 = partial_trace(rho, d, (1, 2))
    rho_f1 = partial_trace(rho, d, (1,))
    return rho_sf1, rho_f1f2, rho_f1, rho


def cmi_nats(rho: Matrix, dims: Tripartition) -> float:
    """Exact von Neumann I(S:F2|F1) in nats."""
    rho_sf1, rho_f1f2, rho_f1, rho_full = marginals(rho, dims)
    return (
        entropy_vn_nats(rho_sf1)
        + entropy_vn_nats(rho_f1f2)
        - entropy_vn_nats(rho_f1)
        - entropy_vn_nats(rho_full)
    )


def i2prime_nats(rho: Matrix, dims: Tripartition) -> float:
    """Surrogate I'2 = H2(SF1)+H2(F1F2)−H2(F1)−H2(SF1F2) in nats.

    Sign-indefinite by design (LMW); never clamped (SEMANTICS §2).
    """
    rho_sf1, rho_f1f2, rho_f1, rho_full = marginals(rho, dims)
    return (
        entropy_renyi2_nats(rho_sf1)
        + entropy_renyi2_nats(rho_f1f2)
        - entropy_renyi2_nats(rho_f1)
        - entropy_renyi2_nats(rho_full)
    )


def sqrtm_psd(rho: Matrix) -> Matrix:
    """ρ^{1/2} for PSD ρ (clip negative numerical noise to 0)."""
    return func_hermitian(rho, np.sqrt, on_support=False)
