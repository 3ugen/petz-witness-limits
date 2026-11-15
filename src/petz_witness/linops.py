"""Tripartite linear algebra: partial traces, embeddings, Hermitian functions.

Index convention (O1): the global space is S ⊗ F1 ⊗ F2 with dimensions
``dims = (d_s, d_1, d_2)``; kron order is S first (most significant).
All density matrices are complex128, trace 1, Hermitian; entropies and
fidelities downstream are in NATS (convention locked in SKELETON header).

Matrix functions of Hermitian operators go through one eigendecomposition
helper so the kernel threshold is applied in exactly one place (O2: one
authority for the support cut).
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "CONF_KERNEL_EPS",
    "Tripartition",
    "func_hermitian",
    "partial_trace",
    "random_density",
    "depolarize",
    "support_projector",
]

Matrix = NDArray[np.complex128]

#: source: Carbon v5 K0 rewrite (reviews/google-carbon-v5/response.md),
#: ratified by owner 2026-11-13 — eigenvalues at or below this are kernel
#: (pseudo-inverted to 0, routed to the CPTP completion), dimensionless.
CONF_KERNEL_EPS = 1e-13


class Tripartition:
    """Dimensions (d_s, d_1, d_2) of S ⊗ F1 ⊗ F2 with derived products."""

    def __init__(self, d_s: int, d_1: int, d_2: int) -> None:
        if d_s < 1 or d_1 < 1 or d_2 < 1:
            raise ValueError(f"dimensions must be >= 1, got {(d_s, d_1, d_2)}")
        self.d_s = d_s
        self.d_1 = d_1
        self.d_2 = d_2

    @property
    def d_total(self) -> int:
        return self.d_s * self.d_1 * self.d_2

    @property
    def d_sf1(self) -> int:
        return self.d_s * self.d_1

    @property
    def d_f1f2(self) -> int:
        return self.d_1 * self.d_2

    def __repr__(self) -> str:
        return f"Tripartition(d_s={self.d_s}, d_1={self.d_1}, d_2={self.d_2})"


def partial_trace(rho: Matrix, dims: tuple[int, ...], keep: tuple[int, ...]) -> Matrix:
    """Partial trace of ρ over the factors NOT in ``keep``.

    ``dims`` are the factor dimensions in kron order; ``keep`` is a sorted
    tuple of factor indices to retain (order of kept factors preserved).
    """
    n = len(dims)
    tensor = rho.reshape(tuple(dims) + tuple(dims))
    traced = [ax for ax in range(n) if ax not in keep]
    # With m remaining row axes the layout is [rows 0..m-1, cols m..2m-1]
    # and row ax pairs with col ax+m; tracing highest ax first keeps the
    # lower axis numbers valid.
    m = n
    for ax in sorted(traced, reverse=True):
        tensor = np.trace(tensor, axis1=ax, axis2=ax + m)
        m -= 1
    d_keep = int(np.prod([dims[i] for i in keep])) if keep else 1
    return np.ascontiguousarray(tensor.reshape(d_keep, d_keep))


def func_hermitian(
    a: Matrix, f: Callable[[NDArray[np.float64]], NDArray[np.float64]], *, on_support: bool
) -> Matrix:
    """f(A) for Hermitian A via eigh; kernel eigenvalues handled in one place.

    ``on_support=True``: eigenvalues ≤ CONF_KERNEL_EPS are mapped to 0
    (pseudo-function on the support — required for ρ^{-1/2}).
    ``on_support=False``: f is applied to eigenvalues clipped at 0
    (numerical noise floor for ρ^{1/2} and spectra).
    """
    w, v = np.linalg.eigh(a)
    w_clipped = np.clip(w, 0.0, None)
    fw = np.zeros_like(w_clipped)
    if on_support:
        mask = w_clipped > CONF_KERNEL_EPS
        fw[mask] = f(w_clipped[mask])
    else:
        fw = f(w_clipped)
    return (v * fw) @ v.conj().T


def support_projector(a: Matrix) -> Matrix:
    """Projector onto the support of Hermitian A (eigenvalues > CONF_KERNEL_EPS)."""
    return func_hermitian(a, lambda w: np.ones_like(w), on_support=True)


def random_density(dim: int, rank: int, rng: np.random.Generator) -> Matrix:
    """Seeded random density matrix of the given rank (Ginibre ρ = GG†/Tr).

    Rank-deficient states (rank < dim) are the K0 anchor class for the
    CPTP completion path; full rank exercises the generic branch.
    """
    if not 1 <= rank <= dim:
        raise ValueError(f"rank must be in [1, {dim}], got {rank}")
    g = rng.standard_normal((dim, rank)) + 1j * rng.standard_normal((dim, rank))
    rho = g @ g.conj().T
    return np.ascontiguousarray((rho / np.trace(rho).real).astype(np.complex128))


def depolarize(rho: Matrix, p: float) -> Matrix:
    """(1−p)·ρ + p·I/d — global white noise (K2 grid reuses this; here anchors)."""
    if not 0.0 <= p <= 1.0:
        raise ValueError(f"p must be in [0, 1], got {p}")
    d = rho.shape[0]
    return (1.0 - p) * rho + p * np.eye(d, dtype=np.complex128) / d
