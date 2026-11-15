"""Petz and twirled-Petz recovery maps R_{F1→SF1}, with CPTP completion.

Sources (O22):
- Unrotated Petz transpose map R^P(Y) = ρ_{SF1}^{1/2}(I_S ⊗ ρ_{F1}^{−1/2}
  Y ρ_{F1}^{−1/2})ρ_{SF1}^{1/2}: Petz 1986/88 as used in BSW 1403.6102
  §8.4 and Fawzi–Renner 1410.0664 (recovery F1 → S F1, N = Tr_S).
- Rotated maps R^{t/2} and the universal averaged ("twirled") map
  R̄ = ∫ β0(t) R^{t/2} dt, β0(t) = (π/2)(cosh(πt)+1)^{−1}: JRSWW
  1509.07127 Theorem 2.1.
- Closed-form spectral kernel ∫ β0(t) e^{iωt} dt = ω/sinh(ω) with the
  PHASE frequency ω[a,x,b,y] = ½ ln(λ_a μ_y / (λ_b μ_x)) (the ½ comes
  from the ρ^{(1±it)/2} powers). Carbon v5 §3 states the same kernel in
  the full-log convention ω' = 2ω as (ω'/2)/sinh(ω'/2); transcribing it
  against the half-log ω here initially produced (ω/2)/sinh(ω/2) — the
  independent quadrature path caught the factor (see test gate), which
  is exactly why two paths exist (T4/O2).
- CPTP completion outside supp(ρ_{F1}): kernel input weight is routed to
  the fixed output ρ_{SF1} (completion freedom per Wilde, Proc. R. Soc. A
  471, 20150338; semantically neutral for actual marginals, whose F1
  support lies inside supp(ρ_{F1}) — tested on rank-deficient anchors).

Consumers: K0 cross-check gate (agreement ≤ 1e-10), K2 oracles. Maps are
produced as Choi matrices of the F1 → SF1 channel in the COMPUTATIONAL
basis of the input factor (O1 convention: J = Σ_{mn} |m⟩⟨n|_in ⊗ M(E_mn),
input factor first), so CP/TP are checked on the object itself and
application to ρ_{F1F2} is one shared contraction.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from petz_witness.entropies import sqrtm_psd
from petz_witness.linops import (
    CONF_KERNEL_EPS,
    Matrix,
    Tripartition,
    func_hermitian,
    partial_trace,
    support_projector,
)

__all__ = [
    "CONF_QUAD_T_MAX",
    "CONF_QUAD_NODES",
    "petz_apply_power",
    "petz_choi",
    "twirled_choi_closed_form",
    "twirled_choi_quadrature",
    "apply_f1_map_choi",
    "assert_cptp",
]

#: ours: β0(t) ≤ (π/2)e^{−π|t|}; at t = 40 the tail weight is ~1e−54 —
#: truncation error far below the 1e−8 cross-check tolerance. Units: t is
#: the dimensionless JRSWW rotation parameter.
CONF_QUAD_T_MAX = 40.0
#: ours: Gauss–Legendre node count. |ω| ≤ ln(1/CONF_KERNEL_EPS) ≈ 30 and
#: |t| ≤ 40 give ≤ ~190 oscillation periods in the window; 4096 nodes are
#: >20 per period. Convergence is then VERIFIED by the agreement gate
#: against the closed form, not assumed (M4).
CONF_QUAD_NODES = 4096


def _beta0(t: NDArray[np.float64]) -> NDArray[np.float64]:
    """JRSWW probability density β0(t) = (π/2)/(cosh(πt)+1)."""
    return (np.pi / 2.0) / (np.cosh(np.pi * t) + 1.0)


def petz_apply_power(
    rho_sf1: Matrix, rho_f1: Matrix, rho_f1f2: Matrix, dims: Tripartition
) -> Matrix:
    """Path 1: R^P(ρ_{F1F2}) by direct matrix powers and embeddings.

    Output on S ⊗ F1 ⊗ F2. Input weight outside supp(ρ_{F1}) goes to
    ρ_{SF1} ⊗ τ_{F2} (τ_{F2} = the F2 remainder of the kernel component),
    matching the Choi path's completion, so the two paths agree on ALL
    inputs, not only on in-support ones.
    """
    d_s, d_2 = dims.d_s, dims.d_2
    inv_sqrt_f1 = func_hermitian(rho_f1, lambda w: 1.0 / np.sqrt(w), on_support=True)
    sqrt_sf1 = sqrtm_psd(rho_sf1)
    pi_f1 = support_projector(rho_f1)
    eye_f2 = np.eye(d_2, dtype=np.complex128)

    inv_f1_emb = np.kron(inv_sqrt_f1, eye_f2)
    y_cond = inv_f1_emb @ rho_f1f2 @ inv_f1_emb
    y_lift = np.kron(np.eye(d_s, dtype=np.complex128), y_cond)
    sqrt_sf1_emb = np.kron(sqrt_sf1, eye_f2)
    out = sqrt_sf1_emb @ y_lift @ sqrt_sf1_emb

    pi_emb = np.kron(pi_f1, eye_f2)
    y_kernel = rho_f1f2 - pi_emb @ rho_f1f2 @ pi_emb
    w_kernel = float(np.real(np.trace(y_kernel)))
    if w_kernel > CONF_KERNEL_EPS:
        tau_f2 = partial_trace(y_kernel, (dims.d_1, d_2), (1,)) / w_kernel
        out = out + np.kron(rho_sf1, tau_f2)
    return out


def _petz_choi_from_weight(
    rho_sf1: Matrix,
    rho_f1: Matrix,
    dims: Tripartition,
    weight: NDArray[np.float64] | None,
) -> Matrix:
    """Choi of the (possibly twirled) completed Petz map F1 → SF1.

    In the eigenbases ρ_{SF1} = Σ λ_a|a⟩⟨a|, ρ_{F1} = Σ μ_x|x⟩⟨x|:
    M(|x⟩⟨y|) = Σ_{a,b} w[a,x,b,y]·amp[a,x]·amp[b,y]
                 · (Σ_s b[a,s,x]·conj(b[b,s,y])-structured) |a⟩⟨b|,
    amp[a,x] = √λ_a/√μ_x, b[a,s,x] = ⟨a|(|s⟩⊗|x⟩). weight=None means the
    unrotated kernel (all ones). The input factor is rotated back to the
    computational basis before the completion term is added, so the Choi
    is consistently computational on both factors.
    """
    d_s, d_1 = dims.d_s, dims.d_1
    d_sf1 = dims.d_sf1
    lam, v_sf1 = np.linalg.eigh(rho_sf1)
    mu, v_f1 = np.linalg.eigh(rho_f1)
    lam = np.clip(lam, 0.0, None)
    mu = np.clip(mu, 0.0, None)
    mu_mask = mu > CONF_KERNEL_EPS

    inv_sqrt_mu = np.where(mu_mask, 1.0 / np.sqrt(np.where(mu_mask, mu, 1.0)), 0.0)
    amp = np.sqrt(lam)[:, None] * inv_sqrt_mu[None, :]  # (a, x)

    # b[a, s, x] = Σ_{f1} conj(v_sf1[(s,f1), a]) · v_f1[f1, x]
    v_sf1_t = v_sf1.reshape(d_s, d_1, d_sf1)
    b = np.einsum("sfa,fx->asx", v_sf1_t.conj(), v_f1)

    # C[a, x, i, s] = v_sf1[i, a] · b[a, s, x]; ca = amp-weighted.
    c = np.einsum("ia,asx->axis", v_sf1, b)
    ca = c * amp[:, :, None, None]

    if weight is None:
        x_blocks = np.einsum("axis,byjs->xiyj", ca, ca.conj())
    else:
        x_blocks = np.einsum("axby,axis,byjs->xiyj", weight, ca, ca.conj())

    # Rotate the input factor to the computational basis:
    # J[m,i,n,j] = Σ_{x,y} conj(v_f1[m,x]) · v_f1[n,y] · X[x,i,y,j]
    j4 = np.einsum("mx,ny,xiyj->minj", v_f1.conj(), v_f1, x_blocks)

    # Completion: kernel directions |x⟩⟨x| (μ_x ≤ eps) → fixed output ρ_{SF1}.
    if not bool(mu_mask.all()):
        rho_out = (v_sf1 * lam) @ v_sf1.conj().T
        for x in range(d_1):
            if not mu_mask[x]:
                vx = v_f1[:, x]
                e_in = np.outer(vx, vx.conj())
                j4 = j4 + np.einsum("mn,ij->minj", e_in, rho_out)

    return np.ascontiguousarray(j4.reshape(d_1 * d_sf1, d_1 * d_sf1))


def petz_choi(rho_sf1: Matrix, rho_f1: Matrix, dims: Tripartition) -> Matrix:
    """Path 2: Choi matrix of the completed unrotated Petz map R^P."""
    return _petz_choi_from_weight(rho_sf1, rho_f1, dims, None)


def _omega(lam: NDArray[np.float64], mu: NDArray[np.float64]) -> NDArray[np.float64]:
    """ω[a,x,b,y] = ½ ln(λ_a μ_y / (λ_b μ_x)) on supports (0 off-support —
    those elements carry amp = 0 and never contribute)."""
    lam_safe = np.where(lam > CONF_KERNEL_EPS, lam, 1.0)
    mu_safe = np.where(mu > CONF_KERNEL_EPS, mu, 1.0)
    ll = np.log(lam_safe)
    lm = np.log(mu_safe)
    return 0.5 * (
        ll[:, None, None, None]
        - ll[None, None, :, None]
        + lm[None, None, None, :]
        - lm[None, :, None, None]
    )


def _twirl_weight_closed(
    lam: NDArray[np.float64], mu: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Closed form: k(ω) = ω/sinh(ω) for the half-log phase ω, k(0) = 1.

    Pinned against scipy.integrate.quad by hand (R-PINNED-BOUND):
    ω = 0.5 → 0.959518, ω = 1.5 → 0.704464, ω = 5 → 0.067383.
    """
    om = _omega(lam, mu)
    out = np.ones_like(om)
    nz = np.abs(om) > 1e-12
    out[nz] = om[nz] / np.sinh(om[nz])
    return out


def _twirl_weight_quadrature(
    lam: NDArray[np.float64], mu: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Independent path: Gauss–Legendre quadrature of ∫ β0(t) cos(ωt) dt
    (the imaginary part vanishes because β0 is even)."""
    nodes, weights = np.polynomial.legendre.leggauss(CONF_QUAD_NODES)
    t = nodes * CONF_QUAD_T_MAX
    w_t = weights * CONF_QUAD_T_MAX * _beta0(t)
    om = _omega(lam, mu)
    return np.tensordot(w_t, np.cos(np.multiply.outer(t, om)), axes=(0, 0))


def twirled_choi_closed_form(rho_sf1: Matrix, rho_f1: Matrix, dims: Tripartition) -> Matrix:
    """Choi of R̄ via the closed-form kernel (primary oracle, SKELETON K0)."""
    lam = np.clip(np.linalg.eigvalsh(rho_sf1), 0.0, None)
    mu = np.clip(np.linalg.eigvalsh(rho_f1), 0.0, None)
    return _petz_choi_from_weight(rho_sf1, rho_f1, dims, _twirl_weight_closed(lam, mu))


def twirled_choi_quadrature(rho_sf1: Matrix, rho_f1: Matrix, dims: Tripartition) -> Matrix:
    """Choi of R̄ via numerical quadrature (independent cross-check path)."""
    lam = np.clip(np.linalg.eigvalsh(rho_sf1), 0.0, None)
    mu = np.clip(np.linalg.eigvalsh(rho_f1), 0.0, None)
    return _petz_choi_from_weight(rho_sf1, rho_f1, dims, _twirl_weight_quadrature(lam, mu))


def apply_f1_map_choi(choi: Matrix, rho_f1f2: Matrix, dims: Tripartition) -> Matrix:
    """(M ⊗ id_{F2})(ρ_{F1F2}) from the Choi of M: F1 → SF1.

    out[(o,a),(p,b)] = Σ_{x,y} J[(x,o),(y,p)] · ρ[(x,a),(y,b)];
    output ordered S ⊗ F1 ⊗ F2 (o = (s,f1) already in that order).
    """
    d_1, d_2 = dims.d_1, dims.d_2
    d_out = dims.d_sf1
    j = choi.reshape(d_1, d_out, d_1, d_out)
    r = rho_f1f2.reshape(d_1, d_2, d_1, d_2)
    out = np.einsum("xoyp,xayb->oapb", j, r)
    return np.ascontiguousarray(out.reshape(d_out * d_2, d_out * d_2))


def assert_cptp(choi: Matrix, dims: Tripartition, *, atol: float) -> None:
    """Raise AssertionError unless the Choi is CP (⪰ −atol) and TP (≤ atol).

    O27 rung: runtime DETECTION on the constructed object (channel-ness).
    What it leaves invisible: map QUALITY — that residue is owed to the
    fidelity oracles and the K0 anchor tests, not to this check.
    """
    w = np.linalg.eigvalsh((choi + choi.conj().T) / 2.0)
    if float(w.min()) < -atol:
        raise AssertionError(f"Choi not CP: min eigenvalue {w.min():.3e}")
    d_1 = dims.d_1
    d_out = dims.d_sf1
    j = choi.reshape(d_1, d_out, d_1, d_out)
    tr_out = np.einsum("xoyo->xy", j)
    dev = float(np.abs(tr_out - np.eye(d_1)).max())
    if dev > atol:
        raise AssertionError(f"Choi not TP: |Tr_out J − I| max {dev:.3e}")
