# petz-witness-limits

**One-sided low-moment witnesses of twirled-Petz recoverability:
noise validity boundaries and shot-budget limits against Pauli state
reconstruction.**

Paper draft + all code, data, and tests. Companion to
[renyi2-darwinism](https://github.com/3ugen/renyi2-darwinism) (the
same methodology applied to Rényi-2 mutual information in quantum
Darwinism); this project asks the conditional-mutual-information /
recoverability question.

## What this is

A sound one-sided certificate U ≥ I(S:F₂|F₁) ≥ −ln F²(ρ, R̄(ρ_{F₁F₂}))
built from second and third Pauli moments of marginals (exact k=2
entropy-gap extremizations, closed-form k=3 ceiling via the
size-biased spectral variance σ² = p₃ − p₂²), its validity map under
noise (9,680 exact-moment evaluations), and a pre-registered negative
finite-shot benchmark: on identical random-Pauli data the witness
loses to projected-MLE state reconstruction by a factor 2–3 at m = 4
and neither method certifies in the sub-tomographic regime. The
headline deliverables are the limits, measured honestly.

- **Paper:** [PAPER.md](PAPER.md) (draft v0.2)
- **Theory with proofs:** [theory/K1_CERTIFICATE.md](theory/K1_CERTIFICATE.md)
- **Gate reports (pre-registered KPI protocol):** [results/](results/)
  — K0 oracles, K1 certificate, K2 validity map, K3 budgets,
  K4 novelty ledger; raw CSVs regenerate every number in the paper.
- **Project log:** [SKELETON.md](SKELETON.md) (gates, statuses,
  correction journal), [SEMANTICS.md](SEMANTICS.md) (scope contract).
- **Adversarial review provenance:** the design and the manuscript
  went through multiple rounds of LLM adversarial review; every
  adopted finding was verified numerically first (analysis notes in
  the parent research log; review prompts under the repository
  history).

## Reproduce

```bash
uv sync --locked
uv run pyrefly check      # 0 errors
uv run pytest -q          # 54 passed
# validity map (~2 min):        uv run python scripts/run_k2.py --fast
# budgets record run (~40 min): uv run python scripts/run_k3.py --m 4 --v4
# near-threshold negatives:     uv run python scripts/run_k3_negatives.py
# figures:                      uv run python scripts/make_figures.py
```

Python ≥ 3.13; numpy/scipy/cvxpy (CLARABEL with automatic SCS fallback
— see the macOS arm64 note in `src/petz_witness/sdp.py`).

## Package map

| Module | Contents |
|---|---|
| `petz_witness.linops` | tripartitions, partial traces, Hermitian calculus |
| `petz_witness.entropies` | CMI, I′₂, moments (nats) |
| `petz_witness.petz` | Petz R^P and twirled R̄, two independent paths each, CPTP completion |
| `petz_witness.certificate` | Lemmas 1–4 machinery: exact δ̄⁽²⁾/δ̲⁽²⁾, closed-form δ̄⁽³⁾, witness U |
| `petz_witness.sampling` | random-Pauli snapshots, exact complete U-statistics, moment-box UCB (v4) |
| `petz_witness.sdp` | CPTP fidelity-of-recovery reference oracle |
| `petz_witness.ensembles` | the eight seeded validity-map ensembles |

MIT license.
