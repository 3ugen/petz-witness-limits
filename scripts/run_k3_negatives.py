"""Run 2 (Carbon v6 §4): near-threshold negative soundness comparison.

Family ρ(λ) = (1−λ)|GHZ₄⟩⟨GHZ₄| + λ·ρ_GHZ-cl,4 with λ tuned by bisection
so the exact twirled-Petz infidelity −ln F_sq(ρ, R̄(ρ_{F1F2})) hits the
targets {0.13, 0.18} — just above τ = 0.10 (true negatives, close to
the decision boundary). Pipelines, all at δ = 0.05:
  A. ours: Û^(3)± v4 UCB            (certify if UCB ≤ τ)
  B. Pauli-MLE plug-in CMI + delete-block jackknife UCB
  C. FoR-SDP on shadow reconstruction + delete-block jackknife UCB
     (reduced arm: fewer reps/blocks — SDP cost; stated in output)
False certification = UCB ≤ τ on a state whose exact infidelity > τ.
A sound 95% one-sided method may false-certify ≤ 5% of runs.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from petz_witness import (
    apply_f1_map_choi,
    cmi_nats,
    marginals,
    neg_log_fidelity_sq,
    twirled_choi_closed_form,
)
from petz_witness.ensembles import _classical_ghz, _pure
from petz_witness.linops import Tripartition
from petz_witness.sampling import (
    PauliDataset,
    baseline_jackknife_ucb,
    sample_dataset,
    shadow_reconstruction,
    witness_ucb_v4,
)

CONF_TAU = 0.10
CONF_TARGETS = (0.13, 0.18)
CONF_REPS = 100
CONF_BUDGETS = (3_000, 10_000, 30_000)
#: ours: SDP arm reduced for cost — reps and blocks below; stated openly.
CONF_SDP_REPS = 30
CONF_SDP_BLOCKS = 8
CONF_SDP_BUDGET = 10_000
CONF_SEED = 20261114

DIMS = Tripartition(2, 4, 2)


def family(lam: float) -> np.ndarray:
    dim = 16
    g = np.zeros(dim, dtype=np.complex128)
    g[0] = g[dim - 1] = 1.0 / np.sqrt(2.0)
    return (1.0 - lam) * _pure(g) + lam * _classical_ghz(4)


def exact_nlf(rho: np.ndarray) -> float:
    rho_sf1, rho_f1f2, rho_f1, _ = marginals(rho, DIMS)
    choi = twirled_choi_closed_form(rho_sf1, rho_f1, DIMS)
    return neg_log_fidelity_sq(rho, apply_f1_map_choi(choi, rho_f1f2, DIMS))


def tune_lambda(target: float) -> tuple[float, np.ndarray]:
    lo, hi = 0.0, 1.0  # nlf decreasing in λ (λ=1 is the classical MC)
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if exact_nlf(family(mid)) > target:
            lo = mid
        else:
            hi = mid
    lam = (lo + hi) / 2.0
    return lam, family(lam)


def sdp_jackknife_ucb(ds: PauliDataset) -> float:
    from scipy.stats import t as student_t

    from petz_witness.sdp import fidelity_of_recovery_sq

    def nlf_of(sub: PauliDataset) -> float:
        rho_hat = shadow_reconstruction(sub, DIMS)
        f = fidelity_of_recovery_sq(rho_hat, DIMS)
        return float(-np.log(max(f, 1e-12)))

    point = nlf_of(ds)
    idx = np.array_split(np.arange(ds.n_snapshots), CONF_SDP_BLOCKS)
    loo = []
    for b in idx:
        keep = np.ones(ds.n_snapshots, dtype=bool)
        keep[b] = False
        loo.append(nlf_of(PauliDataset(ds.codes[keep], ds.n_qubits)))
    loo_arr = np.array(loo)
    k = len(loo_arr)
    pseudo = k * point - (k - 1) * loo_arr
    est = float(pseudo.mean())
    se = float(pseudo.std(ddof=1) / np.sqrt(k))
    tq = float(student_t.ppf(0.95, k - 1))
    return est + tq * se


def main() -> int:
    import csv
    from pathlib import Path

    out_path = Path(__file__).resolve().parent.parent / "results" / "k3_negatives.csv"
    fh = out_path.open("w", newline="")
    writer = csv.writer(fh)
    writer.writerow(["target", "lam", "exact_nlf", "pipeline", "N", "rep", "ucb"])
    rng = np.random.default_rng(CONF_SEED)
    t0 = time.time()
    for target in CONF_TARGETS:
        lam, rho = tune_lambda(target)
        nlf = exact_nlf(rho)
        cmi = cmi_nats(rho, DIMS)
        print(f"\n=== target −lnF = {target}: λ = {lam:.4f}, exact −lnF = {nlf:.4f}, CMI = {cmi:.4f} ===")
        for n_budget in CONF_BUDGETS:
            fc_a = fc_b = 0
            for rep in range(CONF_REPS):
                ds = sample_dataset(rho, DIMS, n_budget, rng)
                _, ucb_a = witness_ucb_v4(ds, DIMS)
                _, ucb_b = baseline_jackknife_ucb(ds, DIMS)
                fc_a += int(ucb_a <= CONF_TAU)
                fc_b += int(ucb_b <= CONF_TAU)
                writer.writerow([target, f"{lam:.6f}", f"{nlf:.6f}", "ours", n_budget, rep, f"{ucb_a:.6f}"])
                writer.writerow([target, f"{lam:.6f}", f"{nlf:.6f}", "mle", n_budget, rep, f"{ucb_b:.6f}"])
                if rep % 10 == 9:
                    fh.flush()
            print(
                f"  N={n_budget:>6}: false-cert ours {fc_a}/{CONF_REPS}, "
                f"MLE-CMI {fc_b}/{CONF_REPS}   ({time.time() - t0:.0f}s)",
                flush=True,
            )
        fc_c = 0
        for rep in range(CONF_SDP_REPS):
            ds = sample_dataset(rho, DIMS, CONF_SDP_BUDGET, rng)
            u = sdp_jackknife_ucb(ds)
            fc_c += int(u <= CONF_TAU)
            writer.writerow([target, f"{lam:.6f}", f"{nlf:.6f}", "sdp", CONF_SDP_BUDGET, rep, f"{u:.6f}"])
            fh.flush()
        print(
            f"  SDP arm (N={CONF_SDP_BUDGET}, {CONF_SDP_REPS} reps, "
            f"{CONF_SDP_BLOCKS} blocks): false-cert FoR-SDP {fc_c}/{CONF_SDP_REPS}"
        )
    fh.close()
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
