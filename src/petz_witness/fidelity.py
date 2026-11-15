"""Fidelity in the LOCKED project convention: squared/Uhlmann F_sq.

Source (O22): F_sq(ρ,σ) = ‖√ρ√σ‖₁² = (Tr √(√ρ σ √ρ))². Convention locked
in SKELETON header after the Carbon v5 audit: Fawzi–Renner's F is ROOT
fidelity, so under F_sq the JRSWW bound reads
I(S:F2|F1) ≥ −ln F_sq(ρ, R̄(ρ_{F1F2})) — the Bell-state anchor test pins
this factor (consumer: every gate comparing −ln F to entropic bounds).
"""

from __future__ import annotations

import numpy as np

from petz_witness.entropies import sqrtm_psd
from petz_witness.linops import Matrix

__all__ = ["fidelity_sq", "neg_log_fidelity_sq"]


def fidelity_sq(rho: Matrix, sigma: Matrix) -> float:
    """F_sq(ρ,σ) = (Tr √(√ρ σ √ρ))², clipped to [0, 1 + numerical noise]."""
    s = sqrtm_psd(rho)
    inner = s @ sigma @ s
    w = np.linalg.eigvalsh((inner + inner.conj().T) / 2.0)
    w = np.clip(w, 0.0, None)
    root = float(np.sqrt(w).sum())
    return root * root


def neg_log_fidelity_sq(rho: Matrix, sigma: Matrix) -> float:
    """−ln F_sq in nats; +inf when the fidelity underflows to 0."""
    f = fidelity_sq(rho, sigma)
    if f <= 0.0:
        return float("inf")
    return float(-np.log(min(f, 1.0)))
