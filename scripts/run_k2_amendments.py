"""Run 3 (Carbon v6 §2): K2 amendments.

(a) p_c(τ) sensitivity curve for τ ∈ [0.02, 0.25], post-processed from
    the existing results/k2_map.csv (no new simulation).
(b) Exact-state U^(3)± recall on ensemble (a) under two additional
    noise channels: 10:1 dephasing-biased Pauli (p_z = 10·p_x = 10·p_y,
    parameter p = total error probability) and amplitude damping (γ).
    Same recall definition and thresholds as K2.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

from petz_witness import (
    apply_f1_map_choi,
    marginals,
    neg_log_fidelity_sq,
    twirled_choi_closed_form,
)
from petz_witness.certificate import witness_u3pm
from petz_witness.ensembles import ensemble_a_stabilizer_perturbed
from petz_witness.entropies import entropy_vn_nats
from petz_witness.linops import Tripartition, partial_trace

CONF_TAU_GRID = (0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.25)
CONF_RECOVERABLE_CUT = 0.02
CONF_MIN_MI = 0.5 * float(np.log(2.0))
CONF_RECALL_GO = 0.20
CONF_SEEDS = 10
CONF_SEED_BASE = 20261113
CONF_NOISE_GRID = (0.0, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2)

_X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
_Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
_Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
_I2 = np.eye(2, dtype=np.complex128)


def _apply_1q_kraus(rho: np.ndarray, kraus: list[np.ndarray], qubit: int, n: int) -> np.ndarray:
    out = np.zeros_like(rho)
    for k in kraus:
        full = np.array([[1.0]], dtype=np.complex128)
        for q in range(n):
            full = np.kron(full, k if q == qubit else _I2)
        out = out + full @ rho @ full.conj().T
    return out


def dephasing_biased(rho: np.ndarray, n: int, p: float) -> np.ndarray:
    """10:1 biased Pauli: p_z = 10 p_x = 10 p_y, total error prob p."""
    pz = p * 10.0 / 12.0
    px = py = p / 12.0
    k = [
        np.sqrt(1.0 - p) * _I2,
        np.sqrt(px) * _X,
        np.sqrt(py) * _Y,
        np.sqrt(pz) * _Z,
    ]
    out = rho
    for q in range(n):
        out = _apply_1q_kraus(out, k, q, n)
    return out


def amplitude_damping(rho: np.ndarray, n: int, gamma: float) -> np.ndarray:
    k0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=np.complex128)
    k1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=np.complex128)
    out = rho
    for q in range(n):
        out = _apply_1q_kraus(out, [k0, k1], q, n)
    return out


def main() -> int:
    base = Path(__file__).resolve().parent.parent / "results"

    # ---- (a) p_c(tau) from the existing map ----
    rows = list(csv.DictReader(open(base / "k2_map.csv")))
    all_p = sorted({float(r["p"]) for r in rows})
    print("p_c(tau) for ensemble (a), witness U^(3)± (from k2_map.csv):")
    for tau in CONF_TAU_GRID:
        pcs = []
        for s in range(CONF_SEEDS):
            cands = []
            for p in all_p:
                pool = [
                    r
                    for r in rows
                    if r["ensemble"] == "a"
                    and float(r["p"]) == p
                    and int(r["seed"]) == s
                    and float(r["mi_sf1"]) >= CONF_MIN_MI
                    and float(r["nlf_twirled"]) <= CONF_RECOVERABLE_CUT
                ]
                if not pool:
                    continue
                rec = sum(1 for r in pool if float(r["u3pm"]) <= tau) / len(pool)
                if rec >= CONF_RECALL_GO:
                    cands.append(p)
            pcs.append(max(cands) if cands else 0.0)
        mean = float(np.mean(pcs))
        spread = (max(pcs) - min(pcs)) / mean if mean > 0 else 0.0
        print(f"  tau={tau:5.2f}: p_c mean={mean:.4g} spread={spread:.1%}")

    # ---- (b) channel comparison, exact states, ensemble (a) ----
    dims_of = {32: Tripartition(2, 4, 4), 16: Tripartition(2, 4, 2)}
    channels = {
        "depolarizing": None,  # ensemble built-in channel handled below
        "dephasing10:1": dephasing_biased,
        "amp-damping": amplitude_damping,
    }
    from petz_witness.ensembles import apply_per_qubit_depolarizing

    print("\nrecall on ensemble (a) by channel (exact states, tau=0.10):")
    for ch_name, ch in channels.items():
        print(f"  --- {ch_name} ---")
        for p in CONF_NOISE_GRID:
            pool_n = 0
            certified = 0
            for seed_i in range(3):  # 3 seeds enough for exact-state recall
                for st in ensemble_a_stabilizer_perturbed(CONF_SEED_BASE + seed_i):
                    dims = st.dims
                    n = int(round(np.log2(dims.d_total)))
                    if ch is None:
                        rho = apply_per_qubit_depolarizing(st.rho, dims, p)
                    else:
                        rho = ch(st.rho, n, p)
                    d3 = (dims.d_s, dims.d_1, dims.d_2)
                    rho_sf1 = partial_trace(rho, d3, (0, 1))
                    rho_s = partial_trace(rho, d3, (0,))
                    rho_f1 = partial_trace(rho, d3, (1,))
                    mi = (
                        entropy_vn_nats(rho_s)
                        + entropy_vn_nats(rho_f1)
                        - entropy_vn_nats(rho_sf1)
                    )
                    if mi < CONF_MIN_MI:
                        continue
                    _, rho_f1f2, rho_f1m, _ = marginals(rho, dims)
                    choi = twirled_choi_closed_form(rho_sf1, rho_f1m, dims)
                    nlf = neg_log_fidelity_sq(
                        rho, apply_f1_map_choi(choi, rho_f1f2, dims)
                    )
                    if nlf > CONF_RECOVERABLE_CUT:
                        continue
                    pool_n += 1
                    if witness_u3pm(rho, dims) <= 0.10:
                        certified += 1
            rec = certified / pool_n if pool_n else float("nan")
            print(f"    p={p:<6}: pool={pool_n:4d} recall={rec:6.1%}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
