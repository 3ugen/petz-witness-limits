"""K2 ensembles (a)–(h), seeded, with the named noise model.

Noise model (O1, named in advance per SKELETON K2): PER-QUBIT
depolarizing with parameter p — each qubit's reduced channel is
ρ_i → (1−p)ρ_i + p·I/2, Kraus {√(1−3p/4)·I, √(p/4)·X, √(p/4)·Y,
√(p/4)·Z}. Grid: CONF_NOISE_GRID (ratified 2026-11-13). The noise is
applied to the WHOLE tripartite state qubit by qubit; "p" in the map
always means this channel, nothing else.

Ensembles (SKELETON §2-K2, a–h). Each generator is a pure function of
its seed; a state is (label, rho, dims). Qubit layout is always
S-block ⊗ F1-block ⊗ F2-block, contiguous (O1).

Sources (O22): central-spin family reuses the p2-darwin-sim model
(H = Σ_k g_k Z_S Z_k, phases exact); HJPW block form per
Hayden–Jozsa–Petz–Winter CMP 246, 359; Gibbs chains H = Σ(Jx XX +
Jz ZZ + h Z) with β-grid — CMI decay regime cited in SKELETON;
graph/cluster chain = CZ-chain on |+⟩^n (its QMC property is VERIFIED
by a test against the exact CMI oracle, not assumed).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from petz_witness.linops import Matrix, Tripartition, random_density

__all__ = [
    "CONF_NOISE_GRID",
    "K2State",
    "apply_per_qubit_depolarizing",
    "ensemble_a_stabilizer_perturbed",
    "ensemble_b_bounded_rank",
    "ensemble_c_central_spin",
    "ensemble_d_lmw_corner",
    "ensemble_e_gibbs_chain",
    "ensemble_f_hjpw",
    "ensemble_g_wghz",
    "ensemble_h_qec",
    "all_ensembles",
]

#: SKELETON K2 noise grid (ratified by owner 2026-11-13), per-qubit
#: depolarizing parameter, dimensionless. 0.0 is the clean reference.
CONF_NOISE_GRID = (0.0, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2)

_X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
_Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
_Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
_I2 = np.eye(2, dtype=np.complex128)
_H = np.array([[1, 1], [1, -1]], dtype=np.complex128) / np.sqrt(2.0)


@dataclass(frozen=True)
class K2State:
    label: str
    rho: Matrix
    dims: Tripartition


def _n_qubits(dims: Tripartition) -> int:
    n = int(round(np.log2(dims.d_total)))
    if 2**n != dims.d_total:
        raise ValueError(f"non-qubit total dimension {dims.d_total}")
    return n


def _apply_1q(rho: Matrix, op: Matrix, qubit: int, n: int) -> Matrix:
    """U ρ U† for a single-qubit operator on the given qubit (kron order)."""
    full = np.array([[1.0]], dtype=np.complex128)
    for q in range(n):
        full = np.kron(full, op if q == qubit else _I2)
    return full @ rho @ full.conj().T


def apply_per_qubit_depolarizing(rho: Matrix, dims: Tripartition, p: float) -> Matrix:
    """The named K2 noise channel, applied to every qubit."""
    if p == 0.0:
        return rho
    n = _n_qubits(dims)
    out = rho
    k0 = float(np.sqrt(1.0 - 3.0 * p / 4.0))
    kp = float(np.sqrt(p / 4.0))
    for q in range(n):
        acc = (k0 * k0) * out
        for pauli in (_X, _Y, _Z):
            acc = acc + (kp * kp) * _apply_1q(out, pauli, q, n)
        out = acc
    return out


def _pure(psi: NDArray[np.complex128]) -> Matrix:
    return np.outer(psi, psi.conj())


def _cz_chain_state(n: int) -> NDArray[np.complex128]:
    """1D cluster state: CZ on neighbours of |+⟩^n."""
    dim = 1 << n
    psi = np.full(dim, 1.0 / np.sqrt(dim), dtype=np.complex128)
    idx = np.arange(dim)
    for q in range(n - 1):
        b1 = (idx >> (n - 1 - q)) & 1
        b2 = (idx >> (n - 2 - q)) & 1
        psi = psi * np.where((b1 & b2) == 1, -1.0, 1.0)
    return psi


def _classical_ghz(n: int) -> Matrix:
    dim = 1 << n
    rho = np.zeros((dim, dim), dtype=np.complex128)
    rho[0, 0] = 0.5
    rho[dim - 1, dim - 1] = 0.5
    return rho


def _random_entangling_rotation(
    rho: Matrix, dims: Tripartition, eps: float, rng: np.random.Generator
) -> Matrix:
    """exp(−i ε G_{q,q+1}) on every adjacent pair (seeded 2-local G).

    Product 1-qubit rotations would be USELESS here: CMI, I'2 and all
    marginal spectra are invariant under local unitaries (observed on
    the first ensemble draft — a product-rotated cluster state kept
    CMI = 0 to 1e-14). Only entangling perturbations move a state off
    the Markov manifold in an invariant way.
    """
    n = _n_qubits(dims)
    out = rho
    for q in range(n - 1):
        g = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))
        g = (g + g.conj().T) / 2.0
        w, v = np.linalg.eigh(g)
        u2 = (v * np.exp(-1j * eps * w)) @ v.conj().T
        full = np.array([[1.0]], dtype=np.complex128)
        pos = 0
        while pos < n:
            if pos == q:
                full = np.kron(full, u2)
                pos += 2
            else:
                full = np.kron(full, _I2)
                pos += 1
        out = full @ out @ full.conj().T
    return out


def ensemble_a_stabilizer_perturbed(seed: int) -> list[K2State]:
    """(a) stabilizer-Markov bases + coherent ε-perturbation."""
    rng = np.random.default_rng(seed)
    out: list[K2State] = []
    bases: list[tuple[str, Matrix, Tripartition]] = [
        ("ghz-cl-4", _classical_ghz(4), Tripartition(2, 4, 2)),
        ("ghz-cl-5", _classical_ghz(5), Tripartition(2, 4, 4)),
        ("cluster-4", _pure(_cz_chain_state(4)), Tripartition(2, 4, 2)),
        ("cluster-5", _pure(_cz_chain_state(5)), Tripartition(2, 4, 4)),
    ]
    for name, rho, dims in bases:
        for eps in (0.0, 0.05, 0.1, 0.2):
            for rep in range(3 if eps > 0.0 else 1):
                r = (
                    rho
                    if eps == 0.0
                    else _random_entangling_rotation(rho, dims, eps, rng)
                )
                out.append(K2State(f"a:{name}:eps{eps}:r{rep}", r, dims))
    return out


def ensemble_b_bounded_rank(seed: int) -> list[K2State]:
    """(b) seeded Ginibre states of bounded rank (soundness fodder)."""
    rng = np.random.default_rng(seed + 1)
    out: list[K2State] = []
    for dims in (Tripartition(2, 4, 2), Tripartition(2, 4, 4)):
        for rank in (1, 2, dims.d_total // 2, dims.d_total):
            for rep in range(3):
                rho = random_density(dims.d_total, rank, rng)
                out.append(K2State(f"b:rank{rank}:r{rep}", rho, dims))
    return out


def ensemble_c_central_spin(seed: int) -> list[K2State]:
    """(c) central-spin pure dephasing (p2 heritage): S + 4 env qubits,
    F1 = env{0,1}, F2 = env{2,3}; |+⟩ start, exact phase evolution."""
    rng = np.random.default_rng(seed + 2)
    out: list[K2State] = []
    n_env = 4
    dims = Tripartition(2, 4, 4)
    n = n_env + 1
    dim = 1 << n
    idx = np.arange(dim)
    for rep in range(3):
        g = rng.uniform(0.3, 1.5, n_env)
        for t in (0.3, 0.8, 1.5):
            psi = np.full(dim, 1.0 / np.sqrt(dim), dtype=np.complex128)
            sigma_s = 1.0 - 2.0 * ((idx >> (n - 1)) & 1)
            energy = np.zeros(dim)
            for k in range(n_env):
                sigma_k = 1.0 - 2.0 * ((idx >> (n - 2 - k)) & 1)
                energy += g[k] * sigma_s * sigma_k
            psi = psi * np.exp(-1j * t * energy)
            out.append(K2State(f"c:g{rep}:t{t}", _pure(psi), dims))
    return out


def ensemble_d_lmw_corner(seed: int) -> list[K2State]:
    """(d) LMW-adversarial corner: classical branching states with
    pointer bias — the I'2 < 0 family (ρ⁻-parametrized)."""
    rng = np.random.default_rng(seed + 3)
    out: list[K2State] = []
    dims = Tripartition(2, 2, 2)
    for w in (0.3, 0.5, 0.7):
        rho = np.zeros((8, 8), dtype=np.complex128)
        # w: S=0, F1=0, F2 uniform; (1−w): S uniform, F1=1, F2=0
        rho[0b000, 0b000] = w / 2.0
        rho[0b001, 0b001] = w / 2.0
        rho[0b010, 0b010] = (1.0 - w) / 2.0
        rho[0b110, 0b110] = (1.0 - w) / 2.0
        out.append(K2State(f"d:rho-minus:w{w}", rho, dims))
    for rep in range(4):
        # random classical 3-bit distributions biased toward branching
        probs = rng.dirichlet(np.full(8, 0.3))
        rho = np.diag(probs.astype(np.complex128))
        out.append(K2State(f"d:classical:r{rep}", rho, dims))
    return out


def _gibbs(n: int, jx: float, jz: float, h: float, beta: float) -> Matrix:
    dim = 1 << n
    ham = np.zeros((dim, dim), dtype=np.complex128)
    for q in range(n - 1):
        for coeff, pauli in ((jx, _X), (jz, _Z)):
            term = np.array([[1.0]], dtype=np.complex128)
            for pos in range(n):
                if pos == q or pos == q + 1:
                    term = np.kron(term, pauli)
                else:
                    term = np.kron(term, _I2)
            ham = ham + coeff * term
    for q in range(n):
        term = np.array([[1.0]], dtype=np.complex128)
        for pos in range(n):
            term = np.kron(term, _Z if pos == q else _I2)
        ham = ham + h * term
    w, v = np.linalg.eigh(ham)
    w = w - w.min()
    boltz = np.exp(-beta * w)
    boltz /= boltz.sum()
    return (v * boltz) @ v.conj().T


def ensemble_e_gibbs_chain(seed: int) -> list[K2State]:
    """(e) 1D Gibbs chains, contiguous S|F1|F2."""
    rng = np.random.default_rng(seed + 4)
    out: list[K2State] = []
    for rep in range(2):
        jx = float(rng.uniform(0.5, 1.2))
        jz = float(rng.uniform(0.2, 0.8))
        h = float(rng.uniform(0.1, 0.6))
        for n, dims in ((4, Tripartition(2, 4, 2)), (5, Tripartition(2, 4, 4))):
            for beta in (0.5, 1.0, 2.0, 4.0):
                rho = _gibbs(n, jx, jz, h, beta)
                out.append(K2State(f"e:n{n}:b{beta}:r{rep}", rho, dims))
    return out


def ensemble_f_hjpw(seed: int) -> list[K2State]:
    """(f) HJPW non-flat exact QMCs: ρ = Σ_j q_j ρ^j_S ⊗ |j⟩⟨j|_{F1L}
    ⊗ σ^j_{F1R F2} — Petz-exact by structure, spectra deliberately
    skewed (the false-refusal family)."""
    rng = np.random.default_rng(seed + 5)
    out: list[K2State] = []
    dims = Tripartition(2, 4, 2)  # F1 = F1L(1 qubit) ⊗ F1R(1 qubit)
    for rep in range(6):
        q0 = float(rng.uniform(0.3, 0.7))
        # near-orthogonal S blocks keep I(S:F1) above the triviality
        # filter 0.5·ln2 (mix ≤ 0.08: binary-entropy loss stays small)
        mix = float(rng.uniform(0.0, 0.08))
        rho_s0 = np.diag([1.0 - mix, mix]).astype(np.complex128)
        rho_s1 = np.diag([mix, 1.0 - mix]).astype(np.complex128)
        sig0 = random_density(4, int(rng.integers(2, 5)), rng)
        sig1 = random_density(4, int(rng.integers(2, 5)), rng)
        e0 = np.diag([1.0, 0.0]).astype(np.complex128)
        e1 = np.diag([0.0, 1.0]).astype(np.complex128)
        rho = q0 * np.kron(rho_s0, np.kron(e0, sig0)) + (1.0 - q0) * np.kron(
            rho_s1, np.kron(e1, sig1)
        )
        out.append(K2State(f"f:hjpw:q{q0:.2f}:r{rep}", rho, dims))
    return out


def ensemble_g_wghz(seed: int) -> list[K2State]:
    """(g) W–GHZ interpolations cosθ|GHZ⟩ + sinθ|W⟩ (orthogonal pair)."""
    out: list[K2State] = []
    for n, dims in ((4, Tripartition(2, 4, 2)), (5, Tripartition(2, 4, 4))):
        dim = 1 << n
        ghz = np.zeros(dim, dtype=np.complex128)
        ghz[0] = ghz[dim - 1] = 1.0 / np.sqrt(2.0)
        w_vec = np.zeros(dim, dtype=np.complex128)
        for q in range(n):
            w_vec[1 << q] = 1.0 / np.sqrt(n)
        for theta in (0.0, 0.3, 0.6, 1.0, np.pi / 2):
            psi = np.cos(theta) * ghz + np.sin(theta) * w_vec
            psi /= np.linalg.norm(psi)
            out.append(K2State(f"g:n{n}:th{theta:.2f}", _pure(psi), dims))
    return out


def ensemble_h_qec(seed: int) -> list[K2State]:
    """(h) decoder diagnostics: reference qubit + 3-qubit repetition
    code; F1 = 2 kept code qubits, F2 = 1 lost; bit-flip and dephasing
    noise on code qubits before the K2 depolarizing sweep."""
    out: list[K2State] = []
    dims = Tripartition(2, 4, 2)
    dim = 16
    psi = np.zeros(dim, dtype=np.complex128)
    psi[0b0000] = 1.0 / np.sqrt(2.0)
    psi[0b1111] = 1.0 / np.sqrt(2.0)
    base = _pure(psi)
    n = 4
    for p_x in (0.0, 0.05, 0.1):
        rho = base
        for q in (1, 2, 3):  # code qubits
            rho = (1.0 - p_x) * rho + p_x * _apply_1q(rho, _X, q, n)
        for p_z in (0.0, 0.2, 1.0):
            r2 = rho
            for q in (1, 2, 3):
                r2 = (1.0 - p_z / 2.0) * r2 + (p_z / 2.0) * _apply_1q(r2, _Z, q, n)
            out.append(K2State(f"h:px{p_x}:pz{p_z}", r2, dims))
    return out


def all_ensembles(seed: int) -> dict[str, list[K2State]]:
    """All eight named ensembles for one seed."""
    return {
        "a": ensemble_a_stabilizer_perturbed(seed),
        "b": ensemble_b_bounded_rank(seed),
        "c": ensemble_c_central_spin(seed),
        "d": ensemble_d_lmw_corner(seed),
        "e": ensemble_e_gibbs_chain(seed),
        "f": ensemble_f_hjpw(seed),
        "g": ensemble_g_wghz(seed),
        "h": ensemble_h_qec(seed),
    }
