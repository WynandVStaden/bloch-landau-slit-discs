"""Exact polish of all local maxima with r>0.99: KD-tree zoom start, then Newton on the
equidistance system with exact nearest points on the exact curves; then verify inscribed."""
import numpy as np, json
from scipy.spatial import cKDTree
from scipy.optimize import minimize_scalar, fsolve
from scipy import ndimage
from dom import R_CLAIM as R, claimed_slits, zslit
S = claimed_slits()
NSAMP = 40001
XS = np.linspace(-1, 0, NSAMP)
SAMP = [zslit(XS, R, al, tip, b) for (c, al, tip, b) in S]
allp = np.concatenate(SAMP)
tree = cKDTree(np.c_[allp.real, allp.imag])
def dist(X, Y):
    d, i = tree.query(np.c_[X.ravel(), Y.ravel()], workers=-1)
    return np.minimum(d, R - np.hypot(X.ravel(), Y.ravel())).reshape(X.shape)
def dist_slit(p, k):
    c, al, tip, b = S[k]
    d = np.abs(SAMP[k] - p); i = d.argmin()
    lo, hi = XS[max(i - 2, 0)], XS[min(i + 2, NSAMP - 1)]
    f = lambda x: abs(zslit(x, R, al, tip, b) - p)
    r = minimize_scalar(f, bounds=(lo, hi), method="bounded", options=dict(xatol=1e-16))
    return min((f(x), x) for x in (r.x, lo, hi))
def pieces(p):
    return [(R - abs(p), None)] + [dist_slit(p, k) for k in range(len(S))]
h = 0.004
g = np.arange(-R, R + h, h)
X, Y = np.meshgrid(g, g)
d = dist(X, Y)
mf = ndimage.maximum_filter(d, size=7)
ii, jj = np.where((d == mf) & (d > 0.99))
starts = []
for i, j in zip(ii, jj):
    cx, cy, hh = X[i, j], Y[i, j], h
    for it in range(30):
        gg = np.linspace(-2 * hh, 2 * hh, 41)
        XX, YY = np.meshgrid(cx + gg, cy + gg)
        dd = dist(XX, YY)
        k = np.unravel_index(dd.argmax(), dd.shape)
        cx, cy, hh = XX[k], YY[k], hh / 4
        if hh < 1e-12: break
    if all(abs(cx + 1j * cy - s) > 1e-5 for s in starts): starts.append(cx + 1j * cy)
res = []
for p in starts:
    ds = pieces(p); dmin = min(d for d, _ in ds)
    ks = [k for k, (dd, x) in enumerate(ds) if dd < dmin + 1e-6]
    if len(ks) < 3:
        continue   # ridge artefact of the grid, not a strict local maximum
    def F(u):
        q = u[0] + 1j * u[1]
        v = [(R - abs(q)) if k == 0 else dist_slit(q, k - 1)[0] for k in ks[:3]]
        return [v[0] - u[2], v[1] - u[2], v[2] - u[2]]
    sol = fsolve(F, [p.real, p.imag, dmin], xtol=1e-15); q = sol[0] + 1j * sol[1]
    ds = pieces(q); rr = min(d for d, _ in ds)
    desc = []; cps = []
    for k, (dd, x) in enumerate(ds):
        if dd < rr + 1e-9:
            if k == 0:
                desc.append("circle"); cps.append(q / abs(q) * R)
            else:
                desc.append(f"{k-1}{S[k-1][0]}:" + ("TIP" if x <= -1 + 1e-12 else f"x={x:.6f}"))
                cps.append(zslit(x, R, *S[k - 1][1:]))
    # strict local max test: centre inside convex hull of contact points (angles span > 180)
    an = np.sort(np.angle(np.array(cps) - q)); gaps = np.diff(np.append(an, an[0] + 2 * np.pi))
    res.append((rr, q, desc, np.degrees(gaps.max()), len(ks)))
res.sort(key=lambda t: -t[0])
for rr, q, desc, gap, m in res:
    print(f"r={rr:.14f} (r-1={rr-1:+.3e}) |c|={abs(q):.10f} arg={np.degrees(np.angle(q))%360:9.5f}  maxgap={gap:6.1f}deg  {desc}")
json.dump([dict(r=rr, x=q.real, y=q.imag, touch=desc) for rr, q, desc, gap, m in res], open("inradius_maxima.json", "w"), indent=1)
print("largest exact inscribed radius: %.15f" % res[0][0])
