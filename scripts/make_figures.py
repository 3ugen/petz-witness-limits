"""Paper figures from recorded CSVs (no recomputation).

Fig. 1: certified recall vs per-qubit depolarizing p on the
stabilizer-Markov ensemble (a), for U^(2), U^(3), U^(3)± (exact
moments; pool = recoverable & non-trivially correlated). Data:
results/k2_map.csv.
Fig. 2: finite-shot certification recall vs N at m = 4 (witness UCB vs
projected-MLE jackknife UCB, 4 positive states x 6 reps). Data:
results/k3_bench_m4_v4.csv.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent / "results"
TAU = 0.10
MI_MIN = 0.5 * math.log(2.0)
REC_CUT = 0.02


def fig1() -> None:
    rows = list(csv.DictReader(open(BASE / "k2_map.csv")))
    ps = sorted({float(r["p"]) for r in rows})
    series = {"u2": [], "u3": [], "u3pm": []}
    for p in ps:
        pool = [
            r
            for r in rows
            if r["ensemble"] == "a"
            and float(r["p"]) == p
            and float(r["mi_sf1"]) >= MI_MIN
            and float(r["nlf_twirled"]) <= REC_CUT
        ]
        for k in series:
            series[k].append(
                sum(1 for r in pool if float(r[k]) <= TAU) / len(pool) if pool else 0.0
            )
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    labels = {"u2": r"$U^{(2)}$ (purities only)", "u3": r"$U^{(3)}$", "u3pm": r"$U^{(3)\pm}$"}
    styles = {"u2": "s:", "u3": "o--", "u3pm": "D-"}
    for k in ("u2", "u3", "u3pm"):
        ax.plot(ps, [100 * v for v in series[k]], styles[k], label=labels[k], ms=4)
    ax.set_xscale("symlog", linthresh=1e-3)
    ax.set_xlabel(r"per-qubit depolarizing $p$")
    ax.set_ylabel("certified recall (%)")
    ax.set_title(r"Exact-moment recall, stabilizer-Markov ensemble ($\tau=0.10$, $m=4$–$5$)")
    ax.axhline(20, color="gray", lw=0.6, ls=":")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(BASE / "fig1_validity.png", dpi=200)
    print("fig1 done")


def fig2() -> None:
    rows = list(csv.DictReader(open(BASE / "k3_bench_m4_v4.csv")))
    ns = sorted({int(r["N"]) for r in rows})
    ours, base = [], []
    for n in ns:
        pool = [r for r in rows if int(r["N"]) == n and r["role"] == "positive"]
        ours.append(100 * sum(1 for r in pool if float(r["ours_ucb"]) <= TAU) / len(pool))
        base.append(100 * sum(1 for r in pool if float(r["base_ucb"]) <= TAU) / len(pool))
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    ax.plot(ns, ours, "D-", label="moment witness $\\hat U$ (95% UCB)", ms=4)
    ax.plot(ns, base, "o--", label="projected-MLE CMI (jackknife UCB)", ms=4)
    ax.set_xscale("log")
    ax.set_xlabel("random-Pauli snapshots $N$")
    ax.set_ylabel("certification recall (%)")
    ax.set_title(r"Finite-shot recall, $m=4$ (4 positives $\times$ 6 reps)")
    ax.axhline(80, color="gray", lw=0.6, ls=":")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(BASE / "fig2_budgets.png", dpi=200)
    print("fig2 done")


if __name__ == "__main__":
    fig1()
    fig2()
