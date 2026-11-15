# Runbook: Run 2 (near-threshold negatives) on the Mac

## What this is

The ONLY missing computation of the p3-petz-witness project (everything
else is done and recorded in results/). Run 2 measures decision-level
false-certification rates of three pipelines on near-threshold
NON-recoverable states — the evidence Carbon v6 demanded before the
paper may print the soundness-delta claim. It was lost to an OOM kill
on the Linux box; the script now writes results incrementally, so a
crash loses at most 10 repetitions.

## Setup (once)

Requires uv (https://docs.astral.sh/uv/) and ~2 GB free RAM steady.

    cd ~/ws/repos/p3-petz-witness
    uv sync --locked          # installs numpy/scipy/cvxpy (+clarabel)

## Verify BEFORE the run (≈2 min; all three must be clean)

    uv run pyrefly check      # expect: 0 errors
    uv run pytest -q          # expect: 54 passed
    uv run python -c "from petz_witness.sdp import fidelity_of_recovery_sq; \
from petz_witness import rho_minus, QUBIT_TRIPARTITION; \
print(fidelity_of_recovery_sq(rho_minus(), QUBIT_TRIPARTITION))"
    # expect: ~0.9999999  (CLARABEL SDP works on this platform)

If pytest or the SDP check fails on macOS (BLAS/solver differences are
the realistic risks), STOP and report the output — do not run Run 2 on
a platform where the gates are red.

## The run (≈1.5–2.5 h single process; safe to Ctrl-C and rerun)

    nohup uv run python scripts/run_k3_negatives.py > run2.log 2>&1 &
    tail -f run2.log          # progress lines every ~100 reps

No pipe through `tail -n` into the main command (that buffering is how
the first attempt's output was lost). Memory stays modest (m = 4,
d = 16); if macOS kills it anyway, the partial CSV is still valid —
report how far it got.

## What "done" looks like

- run2.log ends with: `wrote .../results/k3_negatives.csv`
- For each of two targets (−lnF = 0.13 and 0.18) and each
  N ∈ {3e3, 1e4, 3e4}: a line `false-cert ours X/100, MLE-CMI Y/100`,
  plus one SDP line `false-cert FoR-SDP Z/30`.
- Expected shape (NOT a pass criterion, just orientation): ours ≤ ~5/100
  (the 95% one-sided design level); MLE-CMI and FoR-SDP noticeably
  above 5% at least at the smaller budgets — that asymmetry is the
  paper's claim. If OURS shows >5% substantially, that is a FINDING
  (report it verbatim, do not rerun to make it go away).

## Bring back

results/k3_negatives.csv + run2.log — commit & push, or scp them back;
the analysis and the paper draft consume exactly these two files.
