"""Robustness checks of the BIE result for the claimed domain."""
import numpy as np
from bie import conformal_radius, G0
from dom import R_CLAIM as R, claimed_slits, zslit
sl = [(al, tip, b) for (c, al, tip, b) in claimed_slits()]
for (n, p, power) in [(8, 10, 2), (8, 12, 2), (5, 20, 2), (8, 16, 4), (7, 16, 2)]:
    r, (t, w, Z, mu) = conformal_radius(R, sl, n=n, p=p, power=power)
    print(f"n={n} p={p} power={power}: r={r:.13f}  min mu={mu.min():.3e}")
# non-uniform panels
r, _ = conformal_radius(R, sl, edges=np.array([0, .03, .1, .2, .35, .5, .62, .75, .85, .93, .98, 1.0]))
print("irregular edges: r=%.13f" % r)
# residual: evaluate g at off-slit points close to / away from slits, check harmonic consistency:
r, (t, w, Z, mu) = conformal_radius(R, sl, n=10, p=16)
W = np.tile(w, len(sl))
def g(z):
    return G0(z, 0j, R) + 0 * 0 - np.sum(G0(z, Z, R) * mu * W) if z != 0 else None
# g(z)+log|z| -> log r as z->0
for z in (1e-3, 1e-3j, -1e-4 + 2e-4j):
    val = -np.log(abs(z)) + np.log(R) - np.sum(G0(z, Z, R) * mu * W) + np.log(abs(z))
    print("z=%s: exp(g+log|z|) = %.12f" % (z, np.exp(val)))
# g should be >0 inside D and ~0 approaching slit points between nodes (smooth quadrature valid at distance)
rng = np.random.default_rng(1)
mn = 1e9
for _ in range(20000):
    z = R * np.sqrt(rng.random()) * np.exp(2j * np.pi * rng.random())
    if np.abs(Z - z).min() < 0.05 or abs(z) < 1e-3: continue
    v = -np.log(abs(z)) + np.log(R) - np.sum(G0(z, Z, R) * mu * W)
    mn = min(mn, v)
print("min of g over random interior points (>=0.05 from slits): %.3e (must be >0)" % mn)
print("total charge on slits (harmonic measure of slits from 0): %.12f" % np.sum(mu * W))
