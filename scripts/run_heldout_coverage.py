"""Held-out coverage check (closes the Appendix C item-2 limitation).

The v4 estimator was selected after four revisions on development data
(master seeds 20261113/20261114). This run evaluates UCB coverage of
the exact U on FRESH seeds never used in development, including the
degenerate-marginal case (cluster F1 = I/4) where the normal quantile
is only approximate. Pass criterion (pre-stated): one-sided 95% UCB
covers exact U in >= 93% of runs overall (binomial slack for 350 reps)
and no state family shows coverage below 88%.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from petz_witness.certificate import witness_u3pm
from petz_witness.ensembles import (
    _classical_ghz,
    _cz_chain_state,
    _pure,
    apply_per_qubit_depolarizing,
)
from petz_witness.linops import Tripartition
from petz_witness.sampling import sample_dataset, witness_ucb_v4

CONF_HELDOUT_SEED = 777_001  # fresh; never used in development
CONF_REPS = 25
CONF_BUDGETS = (10_000, 100_000)

DIMS = Tripartition(2, 4, 2)


def states() -> list[tuple[str, np.ndarray]]:
    out = []
    ghz = _classical_ghz(4)
    cl = _pure(_cz_chain_state(4))
    for p in (0.0, 1e-3, 2e-3):
        out.append((f"ghz-cl:p{p}", apply_per_qubit_depolarizing(ghz, DIMS, p)))
    for p in (0.0, 1e-3):
        out.append((f"cluster:p{p}", apply_per_qubit_depolarizing(cl, DIMS, p)))
    dim = 16
    g = np.zeros(dim, dtype=np.complex128)
    g[0] = g[15] = 1 / np.sqrt(2)
    out.append(("ghz-q", _pure(g)))
    w = np.zeros(dim, dtype=np.complex128)
    for q in range(4):
        w[1 << q] = 0.5
    out.append(("w", _pure(w)))
    return out


def main() -> int:
    rng = np.random.default_rng(CONF_HELDOUT_SEED)
    t0 = time.time()
    total = covered = 0
    worst = []
    for name, rho in states():
        exact = witness_u3pm(rho, DIMS)
        st_tot = st_cov = 0
        for n in CONF_BUDGETS:
            for _ in range(CONF_REPS):
                ds = sample_dataset(rho, DIMS, n, rng)
                _, ucb = witness_ucb_v4(ds, DIMS)
                st_tot += 1
                st_cov += int(ucb >= exact - 1e-9)
        total += st_tot
        covered += st_cov
        worst.append((name, st_cov / st_tot))
        print(f"{name:14s}: exact U={exact:.4f} coverage {st_cov}/{st_tot} ({time.time()-t0:.0f}s)", flush=True)
    print(f"\nOVERALL coverage: {covered}/{total} = {covered/total:.1%}")
    fam_min = min(w for _, w in worst)
    verdict = "PASS" if covered / total >= 0.93 and fam_min >= 0.88 else "FAIL"
    print(f"family minimum: {fam_min:.1%}; VERDICT: {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
