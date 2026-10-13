"""Rigorous-style inradius bracket.
Upper bound: branch and bound on squares. d(p) = dist(p, boundary) is 1-Lipschitz, so on a
square with centre c and half-diagonal delta, sup d <= d(c)+delta.  d(c) is evaluated with the
slits replaced by 40001 samples each, which can only OVERESTIMATE d (conservative).
Lower bound / exact values: at each surviving cluster, Newton-solve the equidistance system
with exact nearest-point distances to the exact curves, then verify the disc is inscribed.
"""
import numpy as np, json
from scipy.spatial import cKDTree
from scipy.optimize import minimize_scalar, fsolve
from dom import R_CLAIM as R, claimed_slits, zslit
S = claimed_slits()
NSAMP = 40001
XS = np.linspace(-1, 0, NSAMP)
SAMP = [zslit(XS, R, al, tip, b) for (c, al, tip, b) in S]
allp = np.concatenate(SAMP)
tree = cKDTree(np.c_[allp.real, allp.imag])
def dsamp(cx, cy):
    d, _ = tree.query(np.c_[cx, cy], workers=-1)
    return np.minimum(d, R - np.hypot(cx, cy))
def dist_slit(p, k):
    c, al, tip, b = S[k]
    d = np.abs(SAMP[k] - p); i = d.argmin()
    lo, hi = XS[max(i - 2, 0)], XS[min(i + 2, NSAMP - 1)]
    f = lambda x: abs(zslit(x, R, al, tip, b) - p)
    r = minimize_scalar(f, bounds=(lo, hi), method="bounded", options=dict(xatol=1e-16))
    cands = [(f(x), x) for x in (r.x, lo, hi)]
    return min(cands)
def pieces(p):
    return [(R - abs(p), None)] + [dist_slit(p, k) for k in range(len(S))]

import sys
TARGET = float(sys.argv[1]); STOP = float(sys.argv[2])
side = 0.05
g = np.arange(-R - side, R + side, side) + side / 2
CX, CY = [a.ravel() for a in np.meshgrid(g, g)]
level = 0
while side > STOP:
    delta = side / np.sqrt(2)
    d = dsamp(CX, CY)
    keep = d + delta > TARGET
    CX, CY = CX[keep], CY[keep]
    ub = (d[keep] + delta).max() if keep.any() else None
    print(f"side={side:.3e} survivors={keep.sum():7d}  current upper bound on R(D)={ub}")
    if not keep.any(): break
    if side/2 <= STOP: break
    # subdivide
    q = side / 4
    CX = np.concatenate([CX - q, CX + q, CX - q, CX + q]); CY = np.concatenate([CY - q, CY - q, CY + q, CY + q])
    side /= 2
# cluster survivors
pts = CX + 1j * CY
cl = []
for p in pts:
    for c in cl:
        if abs(c[0] - p) < 1e-3: c.append(p); break
    else: cl.append([p])
print("\nsurviving clusters (each must be explained by an exact local maximum):", len(cl))
res = []
for c in cl:
    p = np.mean(c)
    ds = pieces(p); dmin = min(d for d, _ in ds)
    ks = [k for k, (dd, x) in enumerate(ds) if dd < dmin + 1e-5]
    def F(u):
        q = u[0] + 1j * u[1]
        v = [(R - abs(q)) if k == 0 else dist_slit(q, k - 1)[0] for k in ks[:3]]
        return [v[0] - u[2], v[1] - u[2], v[2] - u[2]]
    if len(ks) >= 3:
        sol = fsolve(F, [p.real, p.imag, dmin], xtol=1e-15); q = sol[0] + 1j * sol[1]
    else:
        q = p
    ds = pieces(q); rr = min(d for d, _ in ds)
    desc = []
    for k, (dd, x) in enumerate(ds):
        if dd < rr + 1e-9:
            desc.append("circle" if k == 0 else f"{k-1}{S[k-1][0]}:" + ("TIP" if x <= -1 + 1e-12 else f"x={x:.6f}"))
    res.append((rr, q, desc, len(c)))
res.sort(key=lambda t: -t[0])
for rr, q, desc, m in res:
    print(f"  exact inscribed radius {rr:.14f} (r-1={rr-1:+.2e}) centre |c|={abs(q):.10f} arg={np.degrees(np.angle(q))%360:9.5f} deg  {desc}")
print("\nLOWER bound R(D) >= %.14f ; UPPER bound R(D) <= %.10f" % (res[0][0], ub))
json.dump([dict(r=rr, x=q.real, y=q.imag, touch=desc) for rr, q, desc, m in res], open("inradius_maxima.json", "w"), indent=1)
