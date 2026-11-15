"""CPTP-optimal fidelity of recovery via SDP (reference oracle).

Sources (O22):
- Definition F(S;F2|F1)_ρ = max_{R CPTP} F(ρ_{SF1F2}, (R⊗id)(ρ_{F1F2})):
  Seshadreesan–Wilde, PRA 92, 042321 (2015), arXiv:1410.1441.
- SDP formulation of ROOT fidelity max ½Tr(X+X†) s.t.
  [[ρ, X], [X†, σ]] ⪰ 0 (Watrous, "Simple semidefinite programs for
  fidelity", arXiv:1207.5726); multiplicativity / SDP treatment of
  fidelity of recovery: Berta–Tomamichel, IEEE TIT 62, 1758 (2016),
  arXiv:1502.07973. NOTE: 1412.4067 is NOT the SDP source (Carbon v5
  citation corrected in the analysis note).

Consumer: K0/K2 reference column "how much is lost to Petz suboptimality
vs moment relaxation". Output is F_sq (project convention): the SDP
optimizes root fidelity, the result is squared once at the end.

Scale note (O7): the Choi variable is (d_1·d_s·d_1)²; intended for
3–6 qubit tripartitions only (SKELETON K0), enforced by a named bound.
"""

from __future__ import annotations

import cvxpy as cp
import numpy as np

from petz_witness.entropies import marginals
from petz_witness.linops import Matrix, Tripartition

__all__ = ["CONF_SDP_MAX_TOTAL_DIM", "fidelity_of_recovery_sq"]

#: ours: largest total Hilbert dimension the SDP oracle accepts; 2^6 = 64
#: keeps the PSD blocks ≤ 128×128 (seconds on a laptop); larger systems
#: are out of the reference oracle's scope (SKELETON K0: 3–6 qubits).
CONF_SDP_MAX_TOTAL_DIM = 64


def fidelity_of_recovery_sq(
    rho: Matrix, dims: Tripartition, *, solver: str = "CLARABEL"
) -> float:
    """F_sq(S;F2|F1) = (max_{R CPTP} ‖√ρ √((R⊗id)ρ_{F1F2})‖₁)².

    Returns the SQUARED fidelity (project convention).
    """
    if dims.d_total > CONF_SDP_MAX_TOTAL_DIM:
        raise ValueError(
            f"total dim {dims.d_total} exceeds CONF_SDP_MAX_TOTAL_DIM="
            f"{CONF_SDP_MAX_TOTAL_DIM} (reference oracle scope)"
        )
    d_1, d_2 = dims.d_1, dims.d_2
    d_out = dims.d_sf1
    d_tot = dims.d_total
    _, rho_f1f2, _, _ = marginals(rho, dims)

    # Choi variable J of R: F1 → SF1, input factor first: (d_1·d_out)².
    j = cp.Variable((d_1 * d_out, d_1 * d_out), hermitian=True)

    # σ = (R ⊗ id_{F2})(ρ_{F1F2}): assembled blockwise,
    # σ = Σ_{x,y} ρ-block2(x,y) ⊗-structured with J-block(x,y):
    # out[(o,a),(p,b)] = Σ_{x,y} J[(x,o),(y,p)] · ρ[(x,a),(y,b)].
    r4 = rho_f1f2.reshape(d_1, d_2, d_1, d_2)
    rows = []
    for o in range(d_out):
        row = []
        for p in range(d_out):
            # block (o,p) in F2 space: Σ_{x,y} J[x*d_out+o, y*d_out+p] · r4[x,:,y,:]
            expr = None
            for x in range(d_1):
                for y in range(d_1):
                    term = j[x * d_out + o, y * d_out + p] * r4[x, :, y, :]
                    expr = term if expr is None else expr + term
            row.append(expr)
        rows.append(row)
    sigma = cp.bmat(rows)  # ordered (S F1) ⊗ F2 = S ⊗ F1 ⊗ F2

    x_var = cp.Variable((d_tot, d_tot), complex=True)
    rho_const = cp.Constant(rho)
    # Tr_out J built by strided slices (cp.partial_trace's stub rejects
    # multi-factor dims tuples): Tr_out[x,y] = Σ_o J[x·d_out+o, y·d_out+o].
    tr_out = j[0::d_out, 0::d_out]
    for o in range(1, d_out):
        tr_out = tr_out + j[o::d_out, o::d_out]
    constraints = [
        cp.bmat([[rho_const, x_var], [cp.conj(x_var).T, sigma]]) >> 0,
        j >> 0,
        tr_out == np.eye(d_1),
    ]
    objective = cp.Maximize(cp.real(cp.trace(x_var)))
    problem = cp.Problem(objective, constraints)
    try:
        problem.solve(solver=solver)
    except cp.error.SolverError:
        # CLARABEL raises SolverError on rank-deficient shadow
        # reconstructions on macOS arm64 (observed 3/3 on Run 2, while
        # solving exact states fine); SCS solves them. Fall back once.
        if solver != "SCS":
            problem.solve(solver="SCS")
        else:
            raise
    value = problem.value
    if problem.status not in ("optimal", "optimal_inaccurate") or value is None:
        raise RuntimeError(f"fidelity-of-recovery SDP failed: status={problem.status}")
    root = float(value)
    return min(root, 1.0) ** 2
