"""K0 anchor states with independently known properties.

Every anchor documents WHAT is pinned and WHERE the expected value comes
from (R-PINNED-BOUND: the oracle value is derived by hand or from a
closed form, never from this code). Order convention: S ⊗ F1 ⊗ F2.

Sources (O22):
- ρ⁻, ρ⁺: Carbon v5 review §4, verified by hand and numerically in
  notes/2026-11-13-carbon-v5-p3-skeleton-analysis.md — classical Markov
  chains (CMI = 0, Petz-recoverable exactly) with I'2 = ln(8/9) < 0 and
  I'2 = ln(10/9) > 0 respectively.
- GHZ / W closed forms: standard (e.g. Nielsen–Chuang); W-state CMI
  derived from two-qubit marginal spectra {2/3, 1/3}.
- HJPW exact QMC structure: Hayden–Jozsa–Petz–Winter, CMP 246, 359
  (2004) — block-diagonal states S–F1L ⊗ F1R–F2 are exactly recoverable.
"""

from __future__ import annotations

import numpy as np

from petz_witness.linops import Matrix, Tripartition

__all__ = [
    "QUBIT_TRIPARTITION",
    "rho_minus",
    "rho_plus",
    "ghz_quantum",
    "ghz_classical",
    "w_state",
    "product_state",
    "bell_pair_trivial_f1",
]

#: One qubit per factor — the anchor geometry (d_s, d_1, d_2) = (2, 2, 2).
QUBIT_TRIPARTITION = Tripartition(2, 2, 2)


def _diag_from_atoms(atoms: dict[tuple[int, int, int], float]) -> Matrix:
    """Diagonal 3-qubit state from classical atoms {(s, f1, f2): prob}."""
    rho = np.zeros((8, 8), dtype=np.complex128)
    for (s, f1, f2), p in atoms.items():
        idx = (s << 2) | (f1 << 1) | f2
        rho[idx, idx] = p
    return rho


def rho_minus() -> Matrix:
    """Classical MC: ½|00⟩⟨00|_{SF1}⊗I/2 + ½(I/2)_S⊗|10⟩⟨10|_{F1F2}.

    Pinned: CMI = 0, F_sq(ρ, R^P(ρ_{F1F2})) = 1, I'2 = ln(8/9) < 0.
    """
    return _diag_from_atoms(
        {(0, 0, 0): 0.25, (0, 0, 1): 0.25, (0, 1, 0): 0.25, (1, 1, 0): 0.25}
    )


def rho_plus() -> Matrix:
    """Classical MC: ½|000⟩⟨000| + ½(I/2)_S⊗|1⟩⟨1|_{F1}⊗(I/2)_{F2}.

    Pinned: CMI = 0, Petz-exact, I'2 = ln(10/9) > 0,
    δ(SF1)+δ(F1F2) = 2(ln(8/3) − H2) … = 0.169925·ln2 nats (analysis note).
    """
    return _diag_from_atoms(
        {
            (0, 0, 0): 0.5,
            (0, 1, 0): 0.125,
            (0, 1, 1): 0.125,
            (1, 1, 0): 0.125,
            (1, 1, 1): 0.125,
        }
    )


def ghz_quantum() -> Matrix:
    """|GHZ₃⟩ = (|000⟩+|111⟩)/√2. Pinned: CMI = ln 2 = I'2 (flat marginals)."""
    psi = np.zeros(8, dtype=np.complex128)
    psi[0] = 1.0 / np.sqrt(2.0)
    psi[7] = 1.0 / np.sqrt(2.0)
    return np.outer(psi, psi.conj())


def ghz_classical() -> Matrix:
    """½(|000⟩⟨000| + |111⟩⟨111|). Pinned: CMI = 0 (classical MC), I'2 = 0."""
    return _diag_from_atoms({(0, 0, 0): 0.5, (1, 1, 1): 0.5})


def w_state() -> Matrix:
    """|W₃⟩ = (|100⟩+|010⟩+|001⟩)/√3.

    Pinned (closed form): single-qubit marginal spectrum {2/3, 1/3};
    two-qubit marginal spectrum {2/3, 1/3}; so
    CMI = H({2/3,1/3}) = ln 3 − (2/3)ln 2 nats.
    """
    psi = np.zeros(8, dtype=np.complex128)
    for idx in (0b100, 0b010, 0b001):
        psi[idx] = 1.0 / np.sqrt(3.0)
    return np.outer(psi, psi.conj())


def product_state() -> Matrix:
    """|0⟩⟨0| ⊗ ρ_mix ⊗ |+⟩⟨+|. Pinned: CMI = 0, I'2 = 0, Petz-exact."""
    rho_s = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=np.complex128)
    rho_f1 = np.array([[0.7, 0.0], [0.0, 0.3]], dtype=np.complex128)
    plus = np.full((2, 2), 0.5, dtype=np.complex128)
    out: Matrix = np.kron(np.kron(rho_s, rho_f1), plus).astype(np.complex128)
    return out


def bell_pair_trivial_f1() -> Matrix:
    """|Φ⁺⟩_{S F2} ⊗ |0⟩⟨0|_{F1} — the fidelity-convention anchor.

    Pinned (hand check in the analysis note): CMI = I(S:F2) = 2 ln 2;
    R̄(ρ_{F1F2}) = |0⟩⟨0|_{F1}-conditioned recovery yields ρ_S ⊗ ρ_{F2}
    marginal product, F_sq(ρ, ρ_S⊗ρ_{F1}⊗ρ_{F2}) = 1/4, so
    −ln F_sq = 2 ln 2 = CMI exactly — a ROOT-fidelity implementation
    would report −2 ln F = 4 ln 2 > CMI and fail this anchor (the M4
    wrong-implementation this test is built to catch).

    Order note: the Bell pair lives on (S, F2); F1 is the middle factor.
    """
    bell = np.zeros(4, dtype=np.complex128)
    bell[0b00] = 1.0 / np.sqrt(2.0)
    bell[0b11] = 1.0 / np.sqrt(2.0)
    rho_sf2 = np.outer(bell, bell.conj())  # order (S, F2)
    # Reorder to S ⊗ F1 ⊗ F2 with F1 = |0⟩⟨0| in the middle.
    rho = np.zeros((8, 8), dtype=np.complex128)
    t = rho_sf2.reshape(2, 2, 2, 2)  # (s, f2, s', f2')
    for s in range(2):
        for f2 in range(2):
            for sp in range(2):
                for f2p in range(2):
                    rho[(s << 2) | (0 << 1) | f2, (sp << 2) | (0 << 1) | f2p] = t[
                        s, f2, sp, f2p
                    ]
    return rho
