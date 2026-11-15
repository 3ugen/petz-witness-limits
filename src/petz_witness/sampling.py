"""K3 sub-tomographic sampling: random-Pauli-basis snapshots and
blocked U-statistic estimators of marginal moments.

Protocol (O1, named): each snapshot = one shot in an independently
random product Pauli basis (uniform over {X,Y,Z}^n) — the "random Pauli
allocation" arm; N_settings = N_total by construction, always ≪ 3^n
coverage in the regime under test. Classical-shadow inverse channel per
qubit: σ = 3·U†|b⟩⟨b|U − I (HKP, Nat. Phys. 16, 1050 — cited for the
protocol; every numeric factor below is COMPUTED from matrices at
import, not transcribed, per the K0 convention lesson).

Estimators: for a marginal M, Tr(ρ̂_i ρ̂_j)|_M = ∏_{q∈M} PAIR[s_i,q, s_j,q]
and Tr(ρ̂_i ρ̂_j ρ̂_k)|_M = Re ∏_{q∈M} TRIPLE[s_i,q, s_j,q, s_k,q], where
s = 2·basis + bit indexes the six single-qubit outcomes. p̂2 = mean over
distinct ordered pairs (full U-statistic, vectorized); p̂3 = mean over a
seeded subsample of distinct triples (incomplete U-statistic — unbiased,
slightly larger variance; O7 bound CONF_TRIPLE_SUBSAMPLE).

Uncertainty: K disjoint blocks of snapshots → one functional estimate
per block → mean + t_{1−δ,K−1}·SE as the one-sided upper confidence
bound (U-statistics per block are iid across blocks). Moments are
clamped into the physical region (p2 ∈ [1/d, 1], p3 ∈ [p2², p2]) BEFORE
the functional — clamping toward the feasible set can only move the
witness toward refusal on garbage data, never toward false
certification of the U-ceiling terms; the I'2 part is clamped
symmetrically and its direction is covered by the block UCB.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from petz_witness.certificate import delta_max3_bound, delta_min2
from petz_witness.linops import Matrix, Tripartition

__all__ = [
    "CONF_BLOCKS_COUNT",
    "CONF_DELTA_CL",
    "CONF_TRIPLE_SUBSAMPLE",
    "PauliDataset",
    "sample_dataset",
    "estimate_p2_blocks",
    "estimate_p3_blocks",
    "witness_ucb_from_dataset",
    "shadow_reconstruction",
]

#: SKELETON K3 (ratified 2026-11-13): one-sided CL level δ = 0.05.
CONF_DELTA_CL = 0.05
#: ours: disjoint blocks for the UCB; 20 keeps t_{0.95,19} = 1.729 while
#: leaving ≥ 50 snapshots per block at the smallest budget (N = 1e3).
CONF_BLOCKS_COUNT = 20
#: ours: triples per block for the incomplete 3rd-order U-statistic —
#: caps the O(N³) count (O7); variance checked against exact moments in
#: tests.
CONF_TRIPLE_SUBSAMPLE = 20_000

_PAULI_EIGVECS: list[NDArray[np.complex128]] = []


def _single_qubit_states() -> list[NDArray[np.complex128]]:
    """The six shadow outcome states |ψ_{basis,bit}⟩, order s = 2b+bit →
    [X+, X−, Y+, Y−, Z+, Z−] as (basis X=0, Y=1, Z=2; bit 0/1)."""
    x = np.array([[1, 1], [1, -1]], dtype=np.complex128) / np.sqrt(2.0)
    y = np.array([[1, 1], [1j, -1j]], dtype=np.complex128) / np.sqrt(2.0)
    z = np.eye(2, dtype=np.complex128)
    out = []
    for u in (x, y, z):
        for bit in (0, 1):
            out.append(u[:, bit].copy())
    return out


def _factor_tables() -> tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """PAIR[s,s'] = Tr(σ_s σ_s'), TRIPLE[s,s',s''] = Tr(σ_s σ_s' σ_s'')
    computed directly from 2×2 matrices (no hand-derived constants)."""
    states = _single_qubit_states()
    sigmas = [3.0 * np.outer(v, v.conj()) - np.eye(2) for v in states]
    pair = np.zeros((6, 6))
    triple = np.zeros((6, 6, 6), dtype=np.complex128)
    for i in range(6):
        for j in range(6):
            pair[i, j] = float(np.real(np.trace(sigmas[i] @ sigmas[j])))
            for k in range(6):
                triple[i, j, k] = np.trace(sigmas[i] @ sigmas[j] @ sigmas[k])
    return pair, triple


_PAIR_TABLE, _TRIPLE_TABLE = _factor_tables()


class PauliDataset:
    """Snapshots as an (N, n) array of outcome codes s = 2·basis + bit."""

    def __init__(self, codes: NDArray[np.int64], n_qubits: int) -> None:
        self.codes = codes
        self.n_qubits = n_qubits

    @property
    def n_snapshots(self) -> int:
        return int(self.codes.shape[0])


def _probability_table(rho: Matrix, n: int) -> NDArray[np.float64]:
    """p[basis_index, outcome] for all 3^n product Pauli bases."""
    states = _single_qubit_states()
    # per-qubit projectors grouped by basis: proj[b][bit]
    n_bases = 3**n
    dim = 1 << n
    table = np.zeros((n_bases, dim))
    # iterate bases; rotate rho per qubit via tensor contraction
    for bi in range(n_bases):
        digits = []
        x = bi
        for _ in range(n):
            digits.append(x % 3)
            x //= 3
        digits.reverse()  # qubit 0 = most significant (kron order)
        # unitary rows: U[bit, :] = ⟨ψ_{basis,bit}|
        t = rho.reshape((2,) * (2 * n))
        for q, b in enumerate(digits):
            u = np.stack(
                [states[2 * b + 0].conj(), states[2 * b + 1].conj()]
            )  # (2 out, 2 in)
            t = np.tensordot(u, t, axes=([1], [q]))
            t = np.moveaxis(t, 0, q)
            uc = u.conj()
            t = np.tensordot(uc, t, axes=([1], [n + q]))
            t = np.moveaxis(t, 0, n + q)
        probs = np.real(
            np.diagonal(t.reshape(dim, dim))
        )
        table[bi] = np.clip(probs, 0.0, None)
        table[bi] /= table[bi].sum()
    return table


def sample_dataset(
    rho: Matrix, dims: Tripartition, n_snapshots: int, rng: np.random.Generator
) -> PauliDataset:
    """Seeded snapshots: uniform random basis per shot, outcome from the
    exact distribution (precomputed 3^n × 2^n table)."""
    n = int(round(np.log2(dims.d_total)))
    table = _probability_table(rho, n)
    n_bases = table.shape[0]
    bases = rng.integers(0, n_bases, n_snapshots)
    outcomes = np.empty(n_snapshots, dtype=np.int64)
    # group by basis for vectorized sampling
    for bi in np.unique(bases):
        mask = bases == bi
        outcomes[mask] = rng.choice(
            table.shape[1], size=int(mask.sum()), p=table[bi]
        )
    # decode into per-qubit codes s = 2*basis_digit + bit
    codes = np.empty((n_snapshots, n), dtype=np.int64)
    for q in range(n):
        digit = (bases // (3 ** (n - 1 - q))) % 3
        bit = (outcomes >> (n - 1 - q)) & 1
        codes[:, q] = 2 * digit + bit
    return PauliDataset(codes, n)


def _pair_products(codes: NDArray[np.int64], qubits: tuple[int, ...]) -> float:
    """Full 2nd-order U-statistic of ∏_{q∈M} PAIR over distinct pairs."""
    n = codes.shape[0]
    prod = np.ones((n, n))
    for q in qubits:
        prod *= _PAIR_TABLE[np.ix_(codes[:, q], codes[:, q])]
    off_diag_sum = float(prod.sum() - np.trace(prod))
    return off_diag_sum / (n * (n - 1))


def _triple_products(
    codes: NDArray[np.int64],
    qubits: tuple[int, ...],
    rng: np.random.Generator,
) -> float:
    """Incomplete 3rd-order U-statistic over a seeded triple subsample."""
    n = codes.shape[0]
    m = min(CONF_TRIPLE_SUBSAMPLE, n * (n - 1) * (n - 2))
    i = rng.integers(0, n, m)
    j = rng.integers(0, n, m)
    k = rng.integers(0, n, m)
    good = (i != j) & (j != k) & (i != k)
    i, j, k = i[good], j[good], k[good]
    prod = np.ones(i.shape[0], dtype=np.complex128)
    for q in qubits:
        prod *= _TRIPLE_TABLE[codes[i, q], codes[j, q], codes[k, q]]
    return float(np.real(prod).mean())


def estimate_p2_blocks(
    ds: PauliDataset, qubits: tuple[int, ...], n_blocks: int = CONF_BLOCKS_COUNT
) -> NDArray[np.float64]:
    """Per-block p̂2 estimates for the marginal on ``qubits``."""
    blocks = np.array_split(np.arange(ds.n_snapshots), n_blocks)
    return np.array(
        [_pair_products(ds.codes[b], qubits) for b in blocks]
    )


def estimate_p3_blocks(
    ds: PauliDataset,
    qubits: tuple[int, ...],
    rng: np.random.Generator,
    n_blocks: int = CONF_BLOCKS_COUNT,
) -> NDArray[np.float64]:
    """Per-block p̂3 estimates (incomplete U-statistic within blocks)."""
    blocks = np.array_split(np.arange(ds.n_snapshots), n_blocks)
    return np.array(
        [_triple_products(ds.codes[b], qubits, rng) for b in blocks]
    )


def _clamp_moments(p2: float, p3: float, d: int) -> tuple[float, float]:
    p2c = min(max(p2, 1.0 / d), 1.0)
    p3c = min(max(p3, p2c * p2c), p2c)
    return p2c, p3c


def witness_ucb_from_dataset(
    ds: PauliDataset,
    dims: Tripartition,
    rng: np.random.Generator,
) -> tuple[float, float]:
    """(point estimate, one-sided 95% UCB) of U^(3)± from snapshots.

    Marginals by qubit position: S = first s-qubits, F1 next, F2 rest.
    """
    from scipy.stats import t as student_t

    n = ds.n_qubits
    n_s = int(round(np.log2(dims.d_s)))
    n_1 = int(round(np.log2(dims.d_1)))
    q_s = tuple(range(n_s))
    q_f1 = tuple(range(n_s, n_s + n_1))
    q_f2 = tuple(range(n_s + n_1, n))
    m_sf1 = q_s + q_f1
    m_f1f2 = q_f1 + q_f2
    m_full = tuple(range(n))

    p2b = {m: estimate_p2_blocks(ds, m) for m in (m_sf1, m_f1f2, q_f1, m_full)}
    p3b = {m: estimate_p3_blocks(ds, m, rng) for m in (m_sf1, m_f1f2)}

    def functional(idx: int) -> float:
        d_sf1, d_f1f2 = dims.d_sf1, dims.d_f1f2
        d_f1, d_tot = dims.d_1, dims.d_total
        p2_sf1, p3_sf1 = _clamp_moments(p2b[m_sf1][idx], p3b[m_sf1][idx], d_sf1)
        p2_f12, p3_f12 = _clamp_moments(p2b[m_f1f2][idx], p3b[m_f1f2][idx], d_f1f2)
        p2_f1, _ = _clamp_moments(p2b[q_f1][idx], 0.0, d_f1)
        p2_full, _ = _clamp_moments(p2b[m_full][idx], 0.0, d_tot)
        i2p = (
            -np.log(p2_sf1) - np.log(p2_f12) + np.log(p2_f1) + np.log(p2_full)
        )
        return float(
            i2p
            + delta_max3_bound(p2_sf1, p3_sf1, d_sf1)
            + delta_max3_bound(p2_f12, p3_f12, d_f1f2)
            - delta_min2(p2_f1, d_f1)
            - delta_min2(p2_full, d_tot)
        )

    k = len(p2b[m_sf1])
    vals = np.array([functional(i) for i in range(k)])
    mean = float(vals.mean())
    se = float(vals.std(ddof=1) / np.sqrt(k))
    tq = float(student_t.ppf(1.0 - CONF_DELTA_CL, k - 1))
    return mean, mean + tq * se


def shadow_reconstruction(ds: PauliDataset, dims: Tripartition) -> Matrix:
    """Mean classical shadow, projected to PSD trace-1 (baseline input)."""
    states = _single_qubit_states()
    sigmas = [3.0 * np.outer(v, v.conj()) - np.eye(2) for v in states]
    n = ds.n_qubits
    dim = dims.d_total
    acc = np.zeros((dim, dim), dtype=np.complex128)
    uniq, counts = np.unique(ds.codes, axis=0, return_counts=True)
    for row, c in zip(uniq, counts):
        term = np.array([[1.0]], dtype=np.complex128)
        for q in range(n):
            term = np.kron(term, sigmas[int(row[q])])
        acc += float(c) * term
    acc /= ds.n_snapshots
    acc = (acc + acc.conj().T) / 2.0
    w, v = np.linalg.eigh(acc)
    w = np.clip(w, 0.0, None)
    s = w.sum()
    if s <= 0.0:
        return np.eye(dim, dtype=np.complex128) / dim
    return (v * (w / s)) @ v.conj().T


# ---------------------------------------------------------------------------
# v2 estimator (2026-11-13, registered after the m=4/m=6 blocked runs):
# block-averaging the NONLINEAR witness functional carries a Jensen bias
# that does not vanish with N (measured: ~1.3 nats at exact 0, m=4,
# N=3e4). v2 estimates MOMENTS on the full dataset (incomplete
# U-statistics with per-snapshot variance components) and propagates a
# joint moment confidence BOX through the functional: upper end for
# terms increasing in a moment, lower end for decreasing, 1D grid over
# the CI where a term is non-monotone (δ ceilings in p2). Conservative
# by construction given the box covers the true moments (Bonferroni
# over 6 moments, two-sided 1 − δ/6 each, normal quantile — U-statistic
# asymptotics, stated caveat for small N).
# ---------------------------------------------------------------------------

#: ours: caps for incomplete U-statistics on the full dataset (O7);
#: variance contributions tracked explicitly, see _moment_with_se.
CONF_PAIR_SUBSAMPLE_FULL = 20_000_000
CONF_TRIPLE_SUBSAMPLE_FULL = 5_000_000
#: ours: chunk size for index subsampling (memory bound ~100 MB).
_CHUNK = 1_000_000


def _moment_with_se(
    codes: "NDArray[np.int64]",
    qubits: tuple[int, ...],
    order: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    """Incomplete U-statistic of the given order with a standard error.

    SE² = Var(h1)·order²/N + Var(kernel)/M  (first term: projection
    variance via per-snapshot mean kernels; second: subsampling noise).
    """
    n = codes.shape[0]
    cap = CONF_PAIR_SUBSAMPLE_FULL if order == 2 else CONF_TRIPLE_SUBSAMPLE_FULL
    total = 0.0
    total_sq = 0.0
    count = 0
    per_snap_sum = np.zeros(n)
    per_snap_cnt = np.zeros(n)
    remaining = cap
    while remaining > 0:
        m = min(_CHUNK, remaining)
        remaining -= m
        idx = [rng.integers(0, n, m) for _ in range(order)]
        good = np.ones(m, dtype=bool)
        for a in range(order):
            for b in range(a + 1, order):
                good &= idx[a] != idx[b]
        idx = [i[good] for i in idx]
        k = idx[0].shape[0]
        if k == 0:
            continue
        if order == 2:
            prod = np.ones(k)
            for q in qubits:
                prod *= _PAIR_TABLE[codes[idx[0], q], codes[idx[1], q]]
        else:
            prodc = np.ones(k, dtype=np.complex128)
            for q in qubits:
                prodc *= _TRIPLE_TABLE[
                    codes[idx[0], q], codes[idx[1], q], codes[idx[2], q]
                ]
            prod = np.real(prodc)
        total += float(prod.sum())
        total_sq += float((prod**2).sum())
        count += k
        np.add.at(per_snap_sum, idx[0], prod)
        np.add.at(per_snap_cnt, idx[0], 1.0)
    mean = total / count
    var_kernel = total_sq / count - mean * mean
    mask = per_snap_cnt >= 3
    h1 = per_snap_sum[mask] / per_snap_cnt[mask]
    var_h1 = float(np.var(h1, ddof=1)) if mask.sum() > 2 else 0.0
    # subtract the kernel-noise inflation of the h1 spread
    var_h1 = max(var_h1 - var_kernel / float(per_snap_cnt[mask].mean()), 0.0)
    se = float(np.sqrt(order * order * var_h1 / n + var_kernel / count))
    return float(mean), se


def witness_ucb_v2(
    ds: PauliDataset,
    dims: Tripartition,
    rng: np.random.Generator,
    delta_cl: float = CONF_DELTA_CL,
) -> tuple[float, float]:
    """(point, one-sided UCB) of U^(3)± via the moment-box propagation."""
    from scipy.stats import norm

    from petz_witness.certificate import delta_max3_bound, delta_min2

    n = ds.n_qubits
    n_s = int(round(np.log2(dims.d_s)))
    n_1 = int(round(np.log2(dims.d_1)))
    q_f1 = tuple(range(n_s, n_s + n_1))
    m_sf1 = tuple(range(n_s)) + q_f1
    m_f1f2 = q_f1 + tuple(range(n_s + n_1, n))
    m_full = tuple(range(n))

    z = float(norm.ppf(1.0 - delta_cl / 12.0))  # 6 moments, two-sided

    est: dict[tuple[tuple[int, ...], int], tuple[float, float]] = {}
    for mq in (m_sf1, m_f1f2, q_f1, m_full):
        est[(mq, 2)] = _moment_with_se(ds.codes, mq, 2, rng)
    for mq in (m_sf1, m_f1f2):
        est[(mq, 3)] = _moment_with_se(ds.codes, mq, 3, rng)

    def box(mq: tuple[int, ...], order: int, d: int) -> tuple[float, float]:
        mean, se = est[(mq, order)]
        lo, hi = mean - z * se, mean + z * se
        lo = min(max(lo, 1.0 / d if order == 2 else 0.0), 1.0)
        hi = min(max(hi, lo), 1.0)
        return lo, hi

    d_sf1, d_f12 = dims.d_sf1, dims.d_f1f2
    d_f1, d_tot = dims.d_1, dims.d_total

    def ceiling_worst(mq: tuple[int, ...], d: int) -> float:
        p2_lo, p2_hi = box(mq, 2, d)
        _, p3_hi = box(mq, 3, d)
        worst = 0.0
        for p2 in np.linspace(p2_lo, p2_hi, 24):
            p2c = min(max(float(p2), 1.0 / d), 1.0)
            p3c = min(max(p3_hi, p2c * p2c), p2c)
            worst = max(worst, delta_max3_bound(p2c, p3c, d))
        return worst

    def min_floor(mq: tuple[int, ...], d: int) -> float:
        p2_lo, p2_hi = box(mq, 2, d)
        return min(
            delta_min2(min(max(float(p2), 1.0 / d), 1.0), d)
            for p2 in np.linspace(p2_lo, p2_hi, 24)
        )

    # I'2 worst case: decreasing in p2(SF1), p2(F1F2); increasing in rest.
    i2_hi = (
        -np.log(box(m_sf1, 2, d_sf1)[0])
        - np.log(box(m_f1f2, 2, d_f12)[0])
        + np.log(box(q_f1, 2, d_f1)[1])
        + np.log(box(m_full, 2, d_tot)[1])
    )
    ucb = float(
        i2_hi
        + ceiling_worst(m_sf1, d_sf1)
        + ceiling_worst(m_f1f2, d_f12)
        - min_floor(q_f1, d_f1)
        - min_floor(m_full, d_tot)
    )

    from petz_witness.certificate import delta_max3_bound as _b  # point est

    def pt_term(mq: tuple[int, ...], d: int) -> float:
        p2, _ = est[(mq, 2)]
        p3, _ = est[(mq, 3)]
        p2c = min(max(p2, 1.0 / d), 1.0)
        p3c = min(max(p3, p2c * p2c), p2c)
        return _b(p2c, p3c, d)

    p2v = {mq: min(max(est[(mq, 2)][0], 1.0 / d), 1.0) for mq, d in
           ((m_sf1, d_sf1), (m_f1f2, d_f12), (q_f1, d_f1), (m_full, d_tot))}
    point = float(
        -np.log(p2v[m_sf1]) - np.log(p2v[m_f1f2])
        + np.log(p2v[q_f1]) + np.log(p2v[m_full])
        + pt_term(m_sf1, d_sf1) + pt_term(m_f1f2, d_f12)
        - delta_min2(p2v[q_f1], d_f1) - delta_min2(p2v[m_full], d_tot)
    )
    return point, ucb


# ---------------------------------------------------------------------------
# v3 (2026-11-13, same session): v2's moment box treats p̂2 and p̂3 as
# independent, which is maximally pessimistic exactly on near-flat
# states, where the projection parts of p̂3 and p̂2 are strongly
# correlated (h1_p3 ≈ p2·h1_p2 on flat spectra, so the delta-method
# variance of σ̂² = p̂3 − p̂2² collapses from 9V to ~V). v3 estimates
# σ² directly with the covariance-aware SE and feeds the ceilings
# (p2-CI, σ²_hi). Bonferroni: 4 two-sided p2 + 2 one-sided σ² → δ/10.
# ---------------------------------------------------------------------------


def _moment_raw(
    codes: "NDArray[np.int64]",
    qubits: tuple[int, ...],
    order: int,
    rng: np.random.Generator,
) -> tuple[float, float, int, "NDArray[np.float64]", "NDArray[np.float64]"]:
    """(mean, kernel variance, count, per-snapshot h1 sums, counts)."""
    n = codes.shape[0]
    cap = CONF_PAIR_SUBSAMPLE_FULL if order == 2 else CONF_TRIPLE_SUBSAMPLE_FULL
    total = total_sq = 0.0
    count = 0
    h1_sum = np.zeros(n)
    h1_cnt = np.zeros(n)
    remaining = cap
    while remaining > 0:
        m = min(_CHUNK, remaining)
        remaining -= m
        idx = [rng.integers(0, n, m) for _ in range(order)]
        good = np.ones(m, dtype=bool)
        for a in range(order):
            for b in range(a + 1, order):
                good &= idx[a] != idx[b]
        idx = [i[good] for i in idx]
        k = idx[0].shape[0]
        if k == 0:
            continue
        if order == 2:
            prod = np.ones(k)
            for q in qubits:
                prod *= _PAIR_TABLE[codes[idx[0], q], codes[idx[1], q]]
        else:
            prodc = np.ones(k, dtype=np.complex128)
            for q in qubits:
                prodc *= _TRIPLE_TABLE[
                    codes[idx[0], q], codes[idx[1], q], codes[idx[2], q]
                ]
            prod = np.real(prodc)
        total += float(prod.sum())
        total_sq += float((prod**2).sum())
        count += k
        np.add.at(h1_sum, idx[0], prod)
        np.add.at(h1_cnt, idx[0], 1.0)
    mean = total / count
    var_kernel = total_sq / count - mean * mean
    return mean, var_kernel, count, h1_sum, h1_cnt


def _sigma2_with_se(
    codes: "NDArray[np.int64]",
    qubits: tuple[int, ...],
    rng: np.random.Generator,
) -> tuple[float, float, float, float]:
    """(p̂2, SE(p̂2), σ̂², SE(σ̂²)) with projection covariance."""
    n = codes.shape[0]
    p2, k2, c2, s2, n2 = _moment_raw(codes, qubits, 2, rng)
    p3, k3, c3, s3, n3 = _moment_raw(codes, qubits, 3, rng)
    mask = (n2 >= 3) & (n3 >= 3)
    h1_2 = s2[mask] / n2[mask]
    h1_3 = s3[mask] / n3[mask]
    v2 = max(float(np.var(h1_2, ddof=1)) - k2 / float(n2[mask].mean()), 0.0)
    v3 = max(float(np.var(h1_3, ddof=1)) - k3 / float(n3[mask].mean()), 0.0)
    c23 = float(np.cov(h1_2, h1_3)[0, 1])  # kernel noises independent
    var_u2 = 4.0 * v2 / n + k2 / c2
    var_u3 = 9.0 * v3 / n + k3 / c3
    cov_u23 = 6.0 * c23 / n
    sigma2 = p3 - p2 * p2
    var_sigma2 = max(
        var_u3 + 4.0 * p2 * p2 * var_u2 - 4.0 * p2 * cov_u23, var_u3 * 1e-4
    )
    return p2, float(np.sqrt(var_u2)), sigma2, float(np.sqrt(var_sigma2))


def witness_ucb_v3(
    ds: PauliDataset,
    dims: Tripartition,
    rng: np.random.Generator,
    delta_cl: float = CONF_DELTA_CL,
) -> tuple[float, float]:
    """(point, one-sided UCB) via covariance-aware σ² propagation."""
    from scipy.stats import norm

    from petz_witness.certificate import delta_max3_bound, delta_min2

    n = ds.n_qubits
    n_s = int(round(np.log2(dims.d_s)))
    n_1 = int(round(np.log2(dims.d_1)))
    q_f1 = tuple(range(n_s, n_s + n_1))
    m_sf1 = tuple(range(n_s)) + q_f1
    m_f1f2 = q_f1 + tuple(range(n_s + n_1, n))
    m_full = tuple(range(n))
    d_of = {m_sf1: dims.d_sf1, m_f1f2: dims.d_f1f2, q_f1: dims.d_1, m_full: dims.d_total}

    z = float(norm.ppf(1.0 - delta_cl / 10.0))

    ceil_est = {mq: _sigma2_with_se(ds.codes, mq, rng) for mq in (m_sf1, m_f1f2)}
    plain_est = {}
    for mq in (q_f1, m_full):
        p2m, kv, cnt, s1, c1 = _moment_raw(ds.codes, mq, 2, rng)
        mask = c1 >= 3
        h1 = s1[mask] / c1[mask]
        v = max(float(np.var(h1, ddof=1)) - kv / float(c1[mask].mean()), 0.0)
        plain_est[mq] = (p2m, float(np.sqrt(4.0 * v / ds.n_snapshots + kv / cnt)))

    def p2_ci(mq: tuple[int, ...]) -> tuple[float, float, float]:
        if mq in ceil_est:
            p2m, se2, _, _ = ceil_est[mq]
        else:
            p2m, se2 = plain_est[mq]
        d = d_of[mq]
        lo = min(max(p2m - z * se2, 1.0 / d), 1.0)
        hi = min(max(p2m + z * se2, lo), 1.0)
        pt = min(max(p2m, 1.0 / d), 1.0)
        return lo, hi, pt

    def ceiling_worst(mq: tuple[int, ...]) -> float:
        d = d_of[mq]
        _, se2, sig2, se_s = ceil_est[mq]
        lo, hi, _ = p2_ci(mq)
        sig_hi = max(sig2 + z * se_s, 0.0)
        worst = 0.0
        for p2 in np.linspace(lo, hi, 16):
            p2c = float(p2)
            p3c = min(max(sig_hi + p2c * p2c, p2c * p2c), p2c)
            worst = max(worst, delta_max3_bound(p2c, p3c, d))
        return worst

    def min_floor(mq: tuple[int, ...]) -> float:
        d = d_of[mq]
        lo, hi, _ = p2_ci(mq)
        return min(
            delta_min2(float(p2), d) for p2 in np.linspace(lo, hi, 16)
        )

    i2_hi = (
        -np.log(p2_ci(m_sf1)[0])
        - np.log(p2_ci(m_f1f2)[0])
        + np.log(p2_ci(q_f1)[1])
        + np.log(p2_ci(m_full)[1])
    )
    ucb = float(
        i2_hi
        + ceiling_worst(m_sf1)
        + ceiling_worst(m_f1f2)
        - min_floor(q_f1)
        - min_floor(m_full)
    )
    # point estimate at the moment point values
    pts = {mq: p2_ci(mq)[2] for mq in (m_sf1, m_f1f2, q_f1, m_full)}

    def pt_ceiling(mq: tuple[int, ...]) -> float:
        d = d_of[mq]
        _, _, sig2, _ = ceil_est[mq]
        p2c = pts[mq]
        p3c = min(max(sig2 + p2c * p2c, p2c * p2c), p2c)
        return delta_max3_bound(p2c, p3c, d)

    point = float(
        -np.log(pts[m_sf1]) - np.log(pts[m_f1f2])
        + np.log(pts[q_f1]) + np.log(pts[m_full])
        + pt_ceiling(m_sf1) + pt_ceiling(m_f1f2)
        - delta_min2(pts[q_f1], d_of[q_f1])
        - delta_min2(pts[m_full], d_of[m_full])
    )
    return point, ucb


def baseline_jackknife_ucb(
    ds: PauliDataset,
    dims: Tripartition,
    n_blocks: int = CONF_BLOCKS_COUNT,
    delta_cl: float = CONF_DELTA_CL,
) -> tuple[float, float]:
    """Full-N shadow/MLE plug-in CMI with delete-block jackknife UCB —
    the strong honest baseline (replaces the per-block average, whose
    reconstruction bias at small blocks dominated)."""
    from scipy.stats import t as student_t

    from petz_witness.entropies import cmi_nats

    full_hat = shadow_reconstruction(ds, dims)
    point = cmi_nats(full_hat, dims)
    idx = np.array_split(np.arange(ds.n_snapshots), n_blocks)
    loo = []
    for b in idx:
        keep = np.ones(ds.n_snapshots, dtype=bool)
        keep[b] = False
        sub = PauliDataset(ds.codes[keep], ds.n_qubits)
        loo.append(cmi_nats(shadow_reconstruction(sub, dims), dims))
    loo_arr = np.array(loo)
    k = len(loo_arr)
    pseudo = k * point - (k - 1) * loo_arr
    est = float(pseudo.mean())
    se = float(pseudo.std(ddof=1) / np.sqrt(k))
    tq = float(student_t.ppf(1.0 - delta_cl, k - 1))
    return est, est + tq * se


# ---------------------------------------------------------------------------
# v4 (2026-11-14, post-Carbon-v6): EXACT complete U-statistics via
# Newton–Girard on accumulated marginal shadows — no subsampling caps
# (the v3 SE(σ̂²)=0.024 was a Monte-Carlo floor from the triple cap,
# verified: kernel var 1724 / 5e6 → 0.019). Formulas verified against
# brute force to 1e-16 (analysis note 2026-11-14, probe B):
#   p̂2 = (N²·Tr ρ̄² − N·5^k) / (N(N−1))
#   p̂3 = (N³·Tr ρ̄³ − 3N·Tr(Q ρ̄) + 2N·7^k) / (N(N−1)(N−2))
# with ρ̄ = mean shadow on M, Q = Σ ρ̂ᵢ², Tr ρ̂ᵢ² = 5^k, Tr ρ̂ᵢ³ = 7^k.
# Exact influence values (Sen):
#   h2_i = (N·Tr(ρ̂ᵢ ρ̄) − 5^k)/(N−1)
#   h3_i = Tr(ρ̂ᵢ[(Nρ̄ − ρ̂ᵢ)² − (Q − ρ̂ᵢ²)])/((N−1)(N−2))
# σ̂² = p̂3 − p̂2² + Var̂(p̂2)  (plug-in bias correction, v6 §3a-3).
# Soundness fixes vs v3: δ_min candidates include every 1/r inside the
# CI (interior zeros — v6 §3a-4, bug confirmed by probe A).
# ---------------------------------------------------------------------------


class _NGMoments:
    """Exact complete U-statistics and influence values for marginal M."""

    def __init__(
        self,
        p2: float,
        p3: float,
        h2: "NDArray[np.float64]",
        h3: "NDArray[np.float64]",
        counts: "NDArray[np.int64]",
        n: int,
    ) -> None:
        self.p2 = p2
        self.p3 = p3
        self.h2 = h2
        self.h3 = h3
        self.counts = counts
        self.n = n


def _marginal_ng(
    codes: "NDArray[np.int64]", qubits: tuple[int, ...]
) -> _NGMoments:
    k = len(qubits)
    dm = 1 << k
    n = codes.shape[0]
    states = _single_qubit_states()
    sigmas = [3.0 * np.outer(v, v.conj()) - np.eye(2) for v in states]

    sub = np.ascontiguousarray(codes[:, list(qubits)])
    uniq, counts = np.unique(sub, axis=0, return_counts=True)

    def build(row: "NDArray[np.int64]") -> "NDArray[np.complex128]":
        t = np.array([[1.0]], dtype=np.complex128)
        for q in range(k):
            t = np.kron(t, sigmas[int(row[q])])
        return t.astype(np.complex128)

    rbar = np.zeros((dm, dm), dtype=np.complex128)
    q_acc = np.zeros((dm, dm), dtype=np.complex128)
    mats_cache: list["NDArray[np.complex128]"] | None = (
        [] if uniq.shape[0] <= 8192 else None
    )
    for row, c in zip(uniq, counts):
        mat = build(row)
        rbar += float(c) * mat
        q_acc += float(c) * (mat @ mat)
        if mats_cache is not None:
            mats_cache.append(mat)
    rbar /= n

    tr2 = float(np.real(np.trace(rbar @ rbar)))
    tr3 = float(np.real(np.trace(rbar @ rbar @ rbar)))
    tr_qr = float(np.real(np.trace(q_acc @ rbar)))
    five_k = 5.0**k
    seven_k = 7.0**k
    p2 = (n * n * tr2 - n * five_k) / (n * (n - 1))
    p3 = (n**3 * tr3 - 3.0 * n * tr_qr + 2.0 * n * seven_k) / (
        n * (n - 1) * (n - 2)
    )

    t_big = n * rbar
    t2_big = t_big @ t_big
    h2 = np.zeros(uniq.shape[0])
    h3 = np.zeros(uniq.shape[0])
    for ui in range(uniq.shape[0]):
        mat = mats_cache[ui] if mats_cache is not None else build(uniq[ui])
        tr_mt2 = float(np.real(np.trace(mat @ t2_big)))
        tr_m2t = float(np.real(np.trace(mat @ mat @ t_big)))
        tr_mq = float(np.real(np.trace(mat @ q_acc)))
        h2[ui] = (n * float(np.real(np.trace(mat @ rbar))) - five_k) / (n - 1)
        h3[ui] = (tr_mt2 - tr_mq - 2.0 * tr_m2t + 2.0 * seven_k) / (
            (n - 1) * (n - 2)
        )
    return _NGMoments(float(p2), float(p3), h2, h3, counts, n)


def _weighted_var_cov(
    a: "NDArray[np.float64]",
    b: "NDArray[np.float64]",
    w: "NDArray[np.int64]",
) -> tuple[float, float, float]:
    """Count-weighted Var(a), Var(b), Cov(a,b)."""
    wt = w / w.sum()
    ma, mb = float((wt * a).sum()), float((wt * b).sum())
    va = float((wt * (a - ma) ** 2).sum())
    vb = float((wt * (b - mb) ** 2).sum())
    cab = float((wt * (a - ma) * (b - mb)).sum())
    return va, vb, cab


def _sigma2_with_se_v4(
    codes: "NDArray[np.int64]", qubits: tuple[int, ...]
) -> tuple[float, float, float, float]:
    """(p̂2, SE(p̂2), σ̂²-corrected, SE(σ̂²)) — exact NG estimators."""
    m = _marginal_ng(codes, qubits)
    n = m.n
    v2, v3, c23 = _weighted_var_cov(m.h2, m.h3, m.counts)
    p2 = m.p2
    p3 = m.p3
    var_u2 = 4.0 * v2 / n
    var_u3 = 9.0 * v3 / n
    cov_u23 = 6.0 * c23 / n
    sigma2 = p3 - p2 * p2 + var_u2  # bias correction E[p̂2²]=p2²+Var
    var_sigma2 = max(
        var_u3 + 4.0 * p2 * p2 * var_u2 - 4.0 * p2 * cov_u23, 0.0
    )
    return p2, float(np.sqrt(var_u2)), sigma2, float(np.sqrt(var_sigma2))


def _p2_candidates(
    lo: float, hi: float, d: int, points: int = 64
) -> "NDArray[np.float64]":
    """Dense grid over [lo, hi] PLUS every flat point 1/r inside — the
    v3 soundness fix (interior zeros of δ_min^(2))."""
    cands = [float(x) for x in np.linspace(lo, hi, points)]
    for r in range(1, d + 1):
        if lo - 1e-15 <= 1.0 / r <= hi + 1e-15:
            cands.append(1.0 / r)
    return np.clip(np.array(cands), 1.0 / d, 1.0)


def witness_ucb_v4(
    ds: PauliDataset,
    dims: Tripartition,
    delta_cl: float = CONF_DELTA_CL,
) -> tuple[float, float]:
    """(point, one-sided UCB) — exact-NG moments, sound candidate sets."""
    from scipy.stats import norm

    from petz_witness.certificate import delta_max3_bound, delta_min2

    n = ds.n_qubits
    n_s = int(round(np.log2(dims.d_s)))
    n_1 = int(round(np.log2(dims.d_1)))
    q_f1 = tuple(range(n_s, n_s + n_1))
    m_sf1 = tuple(range(n_s)) + q_f1
    m_f1f2 = q_f1 + tuple(range(n_s + n_1, n))
    m_full = tuple(range(n))
    d_of = {
        m_sf1: dims.d_sf1,
        m_f1f2: dims.d_f1f2,
        q_f1: dims.d_1,
        m_full: dims.d_total,
    }

    z = float(norm.ppf(1.0 - delta_cl / 10.0))

    ceil_est = {mq: _sigma2_with_se_v4(ds.codes, mq) for mq in (m_sf1, m_f1f2)}
    plain_est = {}
    for mq in (q_f1, m_full):
        mm = _marginal_ng(ds.codes, mq)
        v2, _, _ = _weighted_var_cov(mm.h2, mm.h2, mm.counts)
        plain_est[mq] = (mm.p2, float(np.sqrt(4.0 * v2 / mm.n)))

    def p2_ci(mq: tuple[int, ...]) -> tuple[float, float, float]:
        if mq in ceil_est:
            p2m, se2, _, _ = ceil_est[mq]
        else:
            p2m, se2 = plain_est[mq]
        d = d_of[mq]
        lo = min(max(p2m - z * se2, 1.0 / d), 1.0)
        hi = min(max(p2m + z * se2, lo), 1.0)
        return lo, hi, min(max(p2m, 1.0 / d), 1.0)

    def ceiling_worst(mq: tuple[int, ...]) -> float:
        d = d_of[mq]
        _, _, sig2, se_s = ceil_est[mq]
        lo, hi, _ = p2_ci(mq)
        sig_hi = max(sig2 + z * se_s, 0.0)
        worst = 0.0
        for p2 in _p2_candidates(lo, hi, d):
            p2c = float(p2)
            p3c = min(max(sig_hi + p2c * p2c, p2c * p2c), p2c)
            worst = max(worst, delta_max3_bound(p2c, p3c, d))
        return worst

    def min_floor(mq: tuple[int, ...]) -> float:
        d = d_of[mq]
        lo, hi, _ = p2_ci(mq)
        return min(
            delta_min2(float(p2), d) for p2 in _p2_candidates(lo, hi, d)
        )

    i2_hi = (
        -np.log(p2_ci(m_sf1)[0])
        - np.log(p2_ci(m_f1f2)[0])
        + np.log(p2_ci(q_f1)[1])
        + np.log(p2_ci(m_full)[1])
    )
    ucb = float(
        i2_hi
        + ceiling_worst(m_sf1)
        + ceiling_worst(m_f1f2)
        - min_floor(q_f1)
        - min_floor(m_full)
    )
    pts = {mq: p2_ci(mq)[2] for mq in (m_sf1, m_f1f2, q_f1, m_full)}

    def pt_ceiling(mq: tuple[int, ...]) -> float:
        d = d_of[mq]
        _, _, sig2, _ = ceil_est[mq]
        p2c = pts[mq]
        p3c = min(max(sig2 + p2c * p2c, p2c * p2c), p2c)
        return delta_max3_bound(p2c, p3c, d)

    point = float(
        -np.log(pts[m_sf1])
        - np.log(pts[m_f1f2])
        + np.log(pts[q_f1])
        + np.log(pts[m_full])
        + pt_ceiling(m_sf1)
        + pt_ceiling(m_f1f2)
        - delta_min2(pts[q_f1], d_of[q_f1])
        - delta_min2(pts[m_full], d_of[m_full])
    )
    return point, ucb
