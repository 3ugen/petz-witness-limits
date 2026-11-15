"""K2 validity-map runner (SKELETON §2-K2, thresholds ratified 2026-11-13).

Per (seed, ensemble, state, noise p): records I'2, U^(2), U^(3), exact
CMI, −ln F_sq to twirled R̄ and unrotated R^P, and the filter quantities
I(S:F1), purity moments. Output: results/k2_map.csv (one row per point)
+ stdout summary. Verdict logic (recall/soundness/p_c) lives in the
analysis section at the bottom; thresholds are the ratified [C5] ones.

Run: uv run python scripts/run_k2.py [--seeds 10] [--fast]
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

from petz_witness import (
    apply_f1_map_choi,
    cmi_nats,
    i2prime_nats,
    marginals,
    neg_log_fidelity_sq,
    petz_apply_power,
    twirled_choi_closed_form,
)
from petz_witness.certificate import witness_u2, witness_u3, witness_u3pm
from petz_witness.ensembles import (
    CONF_NOISE_GRID,
    all_ensembles,
    apply_per_qubit_depolarizing,
)
from petz_witness.entropies import entropy_vn_nats
from petz_witness.linops import Tripartition, partial_trace

#: SKELETON K2 ratified thresholds (owner 2026-11-13), nats.
CONF_TAU_NATS = 0.10
CONF_RECOVERABLE_CUT_NATS = 0.02
CONF_TRIVIALITY_MIN_MI_NATS = 0.5 * float(np.log(2.0))
CONF_RECALL_GO = 0.20
CONF_PC_SPREAD_MAX = 0.10
CONF_SEEDS_COUNT = 10
CONF_SEED_BASE = 20261113
#: ours: auxiliary noise points for p_c LOCALIZATION only — the ratified
#: grid stays the gate; extra points add resolution, they cannot weaken
#: the gate (more data on a finer grid).
CONF_AUX_NOISE = (2e-3, 3e-3)


def mutual_info_s_f1(rho: np.ndarray, dims: Tripartition) -> float:
    d = (dims.d_s, dims.d_1, dims.d_2)
    rho_sf1 = partial_trace(rho, d, (0, 1))
    rho_s = partial_trace(rho, d, (0,))
    rho_f1 = partial_trace(rho, d, (1,))
    return (
        entropy_vn_nats(rho_s)
        + entropy_vn_nats(rho_f1)
        - entropy_vn_nats(rho_sf1)
    )


def evaluate_point(
    rho: np.ndarray, dims: Tripartition
) -> dict[str, float]:
    rho_sf1, rho_f1f2, rho_f1, _ = marginals(rho, dims)
    choi = twirled_choi_closed_form(rho_sf1, rho_f1, dims)
    rec_tw = apply_f1_map_choi(choi, rho_f1f2, dims)
    rec_p = petz_apply_power(rho_sf1, rho_f1, rho_f1f2, dims)
    return {
        "i2prime": i2prime_nats(rho, dims),
        "u2": witness_u2(rho, dims),
        "u3": witness_u3(rho, dims),
        "u3pm": witness_u3pm(rho, dims),
        "cmi": cmi_nats(rho, dims),
        "nlf_twirled": neg_log_fidelity_sq(rho, rec_tw),
        "nlf_petz": neg_log_fidelity_sq(rho, rec_p),
        "mi_sf1": mutual_info_s_f1(rho, dims),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=CONF_SEEDS_COUNT)
    parser.add_argument("--fast", action="store_true", help="2 seeds, for smoke runs")
    args = parser.parse_args()
    n_seeds = 2 if args.fast else args.seeds

    out_path = Path(__file__).resolve().parent.parent / "results" / "k2_map.csv"
    out_path.parent.mkdir(exist_ok=True)
    t0 = time.time()
    rows: list[dict[str, object]] = []
    for seed_idx in range(n_seeds):
        seed = CONF_SEED_BASE + seed_idx
        for ens_name, states in all_ensembles(seed).items():
            for st in states:
                for p in tuple(CONF_NOISE_GRID) + CONF_AUX_NOISE:
                    noisy = apply_per_qubit_depolarizing(st.rho, st.dims, p)
                    m = evaluate_point(noisy, st.dims)
                    rows.append(
                        {
                            "seed": seed_idx,
                            "ensemble": ens_name,
                            "label": st.label,
                            "p": p,
                            **{k: f"{v:.10g}" for k, v in m.items()},
                        }
                    )
        print(
            f"seed {seed_idx + 1}/{n_seeds} done, "
            f"{len(rows)} rows, {time.time() - t0:.0f}s",
            flush=True,
        )

    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out_path} ({len(rows)} rows)")

    # ---- analysis: soundness, recall, p_c ----
    data = rows
    float_of = lambda r, k: float(str(r[k]))  # noqa: E731

    viol = [
        r
        for r in data
        if float_of(r, "nlf_twirled")
        > min(float_of(r, "u2"), float_of(r, "u3"), float_of(r, "u3pm")) + 1e-9
    ]
    print(f"\nSOUNDNESS: {len(viol)} violations / {len(data)} points")
    for r in viol[:10]:
        print("  VIOLATION:", r["label"], r["p"], r["u3pm"], r["nlf_twirled"])

    all_p = sorted(set(CONF_NOISE_GRID) | set(CONF_AUX_NOISE))
    print(
        f"\nRECALL per (ensemble, p, witness) — denominators: non-trivial "
        f"(MI >= {CONF_TRIVIALITY_MIN_MI_NATS:.3f}) and recoverable "
        f"(nlf_tw <= {CONF_RECOVERABLE_CUT_NATS}); * marks aux grid points"
    )
    ensembles = sorted({str(r["ensemble"]) for r in data})

    def recall_table(seed_filter: int | None) -> dict[tuple[str, float, str], float]:
        out: dict[tuple[str, float, str], float] = {}
        for ens in ensembles:
            for p in all_p:
                pool = [
                    r
                    for r in data
                    if r["ensemble"] == ens
                    and float(str(r["p"])) == p
                    and (seed_filter is None or r["seed"] == seed_filter)
                    and float_of(r, "mi_sf1") >= CONF_TRIVIALITY_MIN_MI_NATS
                    and float_of(r, "nlf_twirled") <= CONF_RECOVERABLE_CUT_NATS
                ]
                if not pool:
                    continue
                for key in ("u2", "u3", "u3pm"):
                    out[(ens, p, key)] = sum(
                        1 for r in pool if float_of(r, key) <= CONF_TAU_NATS
                    ) / len(pool)
                out[(ens, p, "_pool")] = float(len(pool))
        return out

    table = recall_table(None)
    for ens in ensembles:
        for p in all_p:
            if (ens, p, "_pool") not in table:
                continue
            mark = "*" if p in CONF_AUX_NOISE else " "
            print(
                f"  {ens} p={p:<6}{mark}: pool={int(table[(ens, p, '_pool')]):4d}  "
                f"U2={table[(ens, p, 'u2')]:5.1%}  "
                f"U3={table[(ens, p, 'u3')]:5.1%}  "
                f"U3pm={table[(ens, p, 'u3pm')]:5.1%}"
            )

    # p_c per witness/ensemble: max p with recall >= CONF_RECALL_GO,
    # per-seed spread for the stability clause.
    print(f"\nP_C (max p with recall >= {CONF_RECALL_GO:.0%}) and per-seed spread:")
    for ens in ensembles:
        for key in ("u3", "u3pm"):
            pcs = []
            for s in range(n_seeds):
                t = recall_table(s)
                cands = [
                    p for p in all_p if t.get((ens, p, key), 0.0) >= CONF_RECALL_GO
                ]
                pcs.append(max(cands) if cands else 0.0)
            if max(pcs) == 0.0:
                continue
            mean = float(np.mean(pcs))
            spread = (max(pcs) - min(pcs)) / mean if mean > 0 else 0.0
            print(
                f"  {ens}/{key}: p_c per seed {pcs} -> mean {mean:.4g}, "
                f"rel spread {spread:.1%} (gate <= {CONF_PC_SPREAD_MAX:.0%})"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
