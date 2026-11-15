import numpy as np
from itertools import product

log2 = np.log2

def H1(p):
    p = p[p > 1e-15]
    return -(p * log2(p)).sum()

def H2(p):
    return -log2((p**2).sum())

# --- Counterexample rho^- : classical, diagonal; order (S, F1, F2)
P = {}
# 1/2 |0>_S |0>_F1 (I/2)_F2  -> (0,0,0):1/4, (0,0,1):1/4
P[(0,0,0)] = 0.25; P[(0,0,1)] = 0.25
# 1/2 (I/2)_S |1>_F1 |0>_F2 -> (0,1,0):1/4, (1,1,0):1/4
P[(0,1,0)] = 0.25; P[(1,1,0)] = 0.25

def marg(P, keep):
    out = {}
    for k, v in P.items():
        kk = tuple(k[i] for i in keep)
        out[kk] = out.get(kk, 0) + v
    return np.array(list(out.values()))

def report(P, name):
    full = np.array(list(P.values()))
    sf1  = marg(P, (0,1)); f1f2 = marg(P, (1,2)); f1 = marg(P, (1,))
    I2p = H2(sf1) + H2(f1f2) - H2(f1) - H2(full)
    CMI = H1(sf1) + H1(f1f2) - H1(f1) - H1(full)
    print(f"{name}: I'2 = {I2p:+.6f} bits, CMI = {CMI:+.6f} bits")
    return I2p, CMI

i2m, cmim = report(P, "rho^-")
print("  claim: I'2 = log2(8/9) =", log2(8/9))

# --- Counterexample rho^+
Q = {(0,0,0): 0.5}
for s in (0,1):
    for f2 in (0,1):
        Q[(s,1,f2)] = 0.125
i2p, cmip = report(Q, "rho^+")
print("  claim: I'2 = log2(10/9) =", log2(10/9))
# delta(SF1)+delta(F1F2) for rho^+
sf1 = marg(Q,(0,1)); f1f2 = marg(Q,(1,2))
d = (H1(sf1)-H2(sf1)) + (H1(f1f2)-H2(f1f2))
print(f"  delta(SF1)+delta(F1F2) = {d:.6f} bits (claim ~0.170)")

# --- Counterexample B: rho_SF2 diag((1+e)/4 x2, (1-e)/4 x2), F1 trivial
for eps in (0.1, 0.01, 0.001):
    p = np.array([(1+eps)/4, (1+eps)/4, (1-eps)/4, (1-eps)/4])
    I2 = 1 + 1 - H2(p)            # S2(S)+S2(F2)-S2(SF2), marginals exactly I/2
    I  = 1 + 1 - H1(p)
    Fsq = ((np.sqrt(1+eps) + np.sqrt(1-eps))/2)**2   # fidelity to I/4 (commuting)
    mlnF = -np.log(Fsq)
    I2_nats = I2 * np.log(2); I_nats = I * np.log(2)
    print(f"eps={eps}: I'2/(-lnFsq) = {I2_nats/mlnF:.4f} (claim -> 4), "
          f"I/(-lnFsq) = {I_nats/mlnF:.4f} (claim -> 2)")

# --- marginals of rho_SF2 indeed flat
print("marginal S purity:", 0.5, "-> flat; delta cert slack = 0 (exact)")

# --- depolarized classical GHZ, p=1e-6: CMI and Sum(H0-Hinf)
pdep = 1e-6
# classical 3-bit: 1/2(000)+1/2(111), depolarized: (1-p)*rho + p*I/8 (diagonal)
R = {k: pdep/8 for k in product((0,1),repeat=3)}
R[(0,0,0)] += (1-pdep)/2; R[(1,1,1)] += (1-pdep)/2
i2g, cmig = report(R, "GHZ-cl depol 1e-6")
for keep, name, dim in [((0,1,2),'SF1F2',8),((0,1),'SF1',4),((1,2),'F1F2',4),((1,),'F1',2)]:
    m = marg(R, keep)
    H0 = log2(np.sum(m > 1e-15)); Hinf = -log2(m.max())
    print(f"  {name}: rank={int(2**H0)}, H0-Hinf = {H0-Hinf:.4f} bits")

# --- delta_max^(2)(p2=1/2, d=4) via Harremoes-Topsoe extremal spectrum
d_, p2 = 4, 0.5
x = (1 + np.sqrt((d_-1)*(d_*p2 - 1)))/d_
spec = np.array([x] + [(1-x)/(d_-1)]*(d_-1))
print(f"x = {x:.5f} (claim 0.6830); check p2: {np.sum(spec**2):.6f}")
print(f"delta_max^(2)(1/2,4) = {H1(spec) - H2(spec):.4f} bits (claim 0.4036)")
# brute-force check it is the max of H1 at fixed p2 over 3-level spectra
best = 0
rng = np.random.default_rng(0)
for _ in range(200000):
    v = rng.dirichlet(np.ones(4))
    if abs((v**2).sum() - 0.5) < 2e-3:
        best = max(best, H1(v) - H2(v))
print(f"random search max delta at p2~0.5: {best:.4f} (should be <= bound)")

# --- sigma^2 = p3 - p2^2 >= 0, zero iff flat (spot check)
for spec, nm in [(np.array([0.5,0.5]), "flat r=2"),
                 (np.array([0.5,0.25,0.25]), "non-flat"),
                 (np.ones(3)/3, "flat r=3")]:
    p2_, p3_ = (spec**2).sum(), (spec**3).sum()
    print(f"{nm}: sigma^2 = {p3_ - p2_**2:+.6f}")
