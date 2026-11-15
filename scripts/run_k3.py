"""K3 benchmark: witness-UCB vs same-data shadow/MLE plug-in at m = 6.

Design (SKELETON §2-K3, thresholds ratified; scope notes at bottom of
results/K3_BENCH.md):
- States: 6-qubit tripartition (2, 8, 4), contiguous S|F1|F2.
  Positives = truly recoverable (−lnF_sq(ρ, R̄) ≤ 0.02) AND exactly
  certifiable (U^(3)± ≤ τ = 0.10); negatives = not recoverable
  (−lnF > τ). Borderline states are excluded from recall denominators
  but count for soundness.
- Budgets N ∈ {1e3, 3e3, 1e4, 3e4} snapshots (= settings; random Pauli
  allocation, N_settings ≪ 3^6 = 729 only in the sense of coverage per
  marginal moment — every snapshot is its own setting).
- Ours: Û^(3)±_{95%} block-UCB ≤ τ certifies.
- Baseline (same snapshots): per-block shadow reconstruction (PSD
  projection) → plug-in CMI(ρ̂_block) → block-UCB ≤ τ certifies.
- Metrics: recall on positives per budget; empirical UCB coverage of
  the exact value; decision-level false certifications of negatives;
  GO ratio = N*_baseline / N*_ours at 80% recall.

Run: uv run python scripts/run_k3.py [--fast]
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

from petz_witness import (
    apply_f1_map_choi,
    cmi_nats,
    marginals,
    neg_log_fidelity_sq,
    twirled_choi_closed_form,
)
from petz_witness.certificate import witness_u3pm
from petz_witness.ensembles import (
    _classical_ghz,
    _cz_chain_state,
    _pure,
    apply_per_qubit_depolarizing,
)
from petz_witness.linops import Tripartition
from petz_witness.sampling import (
    CONF_BLOCKS_COUNT,
    CONF_DELTA_CL,
    sample_dataset,
    shadow_reconstruction,
    witness_ucb_from_dataset,
    PauliDataset,
)

#: ratified [C5]
CONF_TAU_NATS = 0.10
CONF_RECOVERABLE_CUT = 0.02
CONF_RECALL_TARGET = 0.80
CONF_GO_RATIO = 1.5
#: ours: benchmark repetition count per (state, budget) — 8 reps bound
#: the recall resolution to 1/8; chosen for wall-clock, stated openly.
CONF_REPS = 8
CONF_BUDGETS = (1_000, 3_000, 10_000, 30_000)
CONF_SEED = 20261113

DIMS = Tripartition(2, 8, 4)


def build_states() -> list[tuple[str, np.ndarray]]:
    n = 6 if DIMS.d_total == 64 else 4
    out = []
    ghz = _classical_ghz(n)
    cluster = _pure(_cz_chain_state(n))
    noise_levels = (0.0, 1e-3, 2e-3) if n == 6 else (0.0, 1e-3)
    for p in noise_levels:
        out.append((f"ghz-cl-{n}:p{p}", apply_per_qubit_depolarizing(ghz, DIMS, p)))
        out.append((f"cluster-{n}:p{p}", apply_per_qubit_depolarizing(cluster, DIMS, p)))
    # negatives: quantum GHZ (CMI = ln 2, not recoverable), W state
    dim = 1 << n
    g = np.zeros(dim, dtype=np.complex128)
    g[0] = g[dim - 1] = 1.0 / np.sqrt(2.0)
    out.append((f"ghz-q-{n}", _pure(g)))
    w = np.zeros(dim, dtype=np.complex128)
    for q in range(n):
        w[1 << q] = 1.0 / np.sqrt(n)
    out.append((f"w-{n}", _pure(w)))
    return out


def exact_quantities(rho: np.ndarray) -> dict[str, float]:
    rho_sf1, rho_f1f2, rho_f1, _ = marginals(rho, DIMS)
    choi = twirled_choi_closed_form(rho_sf1, rho_f1, DIMS)
    rec = apply_f1_map_choi(choi, rho_f1f2, DIMS)
    return {
        "nlf": neg_log_fidelity_sq(rho, rec),
        "u3pm": witness_u3pm(rho, DIMS),
        "cmi": cmi_nats(rho, DIMS),
    }


def baseline_block_ucb(ds: PauliDataset, n_blocks: int = CONF_BLOCKS_COUNT) -> tuple[float, float]:
    """Same-data baseline: per-block shadow/MLE plug-in CMI → block UCB."""
    idx = np.array_split(np.arange(ds.n_snapshots), n_blocks)
    vals = []
    for b in idx:
        sub = PauliDataset(ds.codes[b], ds.n_qubits)
        rho_hat = shadow_reconstruction(sub, DIMS)
        vals.append(cmi_nats(rho_hat, DIMS))
    arr = np.array(vals)
    mean = float(arr.mean())
    se = float(arr.std(ddof=1) / np.sqrt(len(arr)))
    tq = float(student_t.ppf(1.0 - CONF_DELTA_CL, len(arr) - 1))
    return mean, mean + tq * se


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fast", action="store_true")
    parser.add_argument("--m", type=int, default=6, choices=(4, 6))
    parser.add_argument("--v3", action="store_true", help="v3 estimators, record run")
    parser.add_argument("--v4", action="store_true", help="v4 exact-NG estimators (post-Carbon-v6)")
    args = parser.parse_args()
    global DIMS
    if args.m == 4:
        DIMS = Tripartition(2, 4, 2)
    reps = 3 if args.fast else (6 if (args.v3 or args.v4) else CONF_REPS)
    if args.fast:
        budgets = CONF_BUDGETS[:2]
    elif args.v4:
        budgets = (
            (10_000, 30_000, 100_000, 300_000, 1_000_000)
            if args.m == 4
            else (300, 3_000, 30_000, 300_000)
        )
        if args.m == 6:
            reps = 3
    elif args.v3:
        budgets = (30_000, 100_000, 300_000)
    else:
        budgets = CONF_BUDGETS

    states = build_states()
    print("exact classification:")
    classified = []
    for name, rho in states:
        q = exact_quantities(rho)
        role = (
            "positive"
            if q["nlf"] <= CONF_RECOVERABLE_CUT and q["u3pm"] <= CONF_TAU_NATS
            else "negative"
            if q["nlf"] > CONF_TAU_NATS
            else "borderline"
        )
        classified.append((name, rho, q, role))
        print(
            f"  {name:<18} nlf={q['nlf']:.4f} U3pm={q['u3pm']:.4f} "
            f"cmi={q['cmi']:.4f} -> {role}"
        )

    rng = np.random.default_rng(CONF_SEED)
    rows = []
    t0 = time.time()
    for name, rho, q, role in classified:
        for n_budget in budgets:
            for rep in range(reps):
                ds = sample_dataset(rho, DIMS, n_budget, rng)
                if args.v4:
                    from petz_witness.sampling import (
                        baseline_jackknife_ucb,
                        witness_ucb_v4,
                    )

                    ours_pt, ours_ucb = witness_ucb_v4(ds, DIMS)
                    base_pt, base_ucb = baseline_jackknife_ucb(ds, DIMS)
                elif args.v3:
                    from petz_witness.sampling import (
                        baseline_jackknife_ucb,
                        witness_ucb_v3,
                    )

                    ours_pt, ours_ucb = witness_ucb_v3(ds, DIMS, rng)
                    base_pt, base_ucb = baseline_jackknife_ucb(ds, DIMS)
                else:
                    ours_pt, ours_ucb = witness_ucb_from_dataset(ds, DIMS, rng)
                    base_pt, base_ucb = baseline_block_ucb(ds)
                rows.append(
                    {
                        "state": name,
                        "role": role,
                        "N": n_budget,
                        "rep": rep,
                        "ours_pt": ours_pt,
                        "ours_ucb": ours_ucb,
                        "base_pt": base_pt,
                        "base_ucb": base_ucb,
                        "exact_u3pm": q["u3pm"],
                        "exact_cmi": q["cmi"],
                        "exact_nlf": q["nlf"],
                    }
                )
            print(
                f"{name} N={n_budget}: done ({time.time() - t0:.0f}s)", flush=True
            )

    fname = (
        f"k3_bench_m{6 if DIMS.d_total == 64 else 4}"
        + ("_v4" if args.v4 else "_v3" if args.v3 else "")
        + ".csv"
    )
    out_path = Path(__file__).resolve().parent.parent / "results" / fname
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out_path} ({len(rows)} rows)")

    # ---- analysis ----
    def recall(method_key: str, n_budget: int) -> float:
        pool = [r for r in rows if r["role"] == "positive" and r["N"] == n_budget]
        return sum(1 for r in pool if r[method_key] <= CONF_TAU_NATS) / len(pool)

    print("\nRECALL on positives (certify = UCB <= tau):")
    n_star = {"ours_ucb": None, "base_ucb": None}
    for n_budget in budgets:
        ro = recall("ours_ucb", n_budget)
        rb = recall("base_ucb", n_budget)
        for key, val in (("ours_ucb", ro), ("base_ucb", rb)):
            if val >= CONF_RECALL_TARGET and n_star[key] is None:
                n_star[key] = n_budget
        print(f"  N={n_budget:>6}: ours={ro:5.1%}  baseline={rb:5.1%}")

    cov_o = sum(
        1 for r in rows if r["ours_ucb"] >= r["exact_u3pm"] - 1e-9
    ) / len(rows)
    cov_b = sum(1 for r in rows if r["base_ucb"] >= r["exact_cmi"] - 1e-9) / len(
        rows
    )
    fc_o = sum(
        1
        for r in rows
        if r["role"] == "negative" and r["ours_ucb"] <= CONF_TAU_NATS
    )
    fc_b = sum(
        1
        for r in rows
        if r["role"] == "negative" and r["base_ucb"] <= CONF_TAU_NATS
    )
    neg_n = sum(1 for r in rows if r["role"] == "negative")
    print(f"\nUCB coverage of exact value: ours {cov_o:.1%}, baseline {cov_b:.1%}")
    print(f"False certifications of negatives: ours {fc_o}/{neg_n}, baseline {fc_b}/{neg_n}")
    print(f"\nN* (80% recall): ours={n_star['ours_ucb']}, baseline={n_star['base_ucb']}")
    if n_star["ours_ucb"] and n_star["base_ucb"]:
        print(f"GO ratio (base/ours) = {n_star['base_ucb'] / n_star['ours_ucb']:.2f} (gate >= {CONF_GO_RATIO})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
