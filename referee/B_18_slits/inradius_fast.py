"""Fast second inradius pass: fine sampling (40001 pts/slit, spacing<=1e-4 => distance
overestimate <= ~2e-9 for d~1), global grid h=0.004, then local nested-grid zoom using KD-tree."""
import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage
from dom import R_CLAIM as R, claimed_slits, zslit
S = claimed_slits()
XS = np.linspace(-1, 0, 40001)
allp = np.concatenate([zslit(XS, R, al, tip, b) for (c, al, tip, b) in S])
lab = np.repeat(np.arange(len(S)), len(XS))
tree = cKDTree(np.c_[allp.real, allp.imag])
def dist(X, Y):
    d, i = tree.query(np.c_[X.ravel(), Y.ravel()], workers=-1)
    return np.minimum(d, R - np.hypot(X.ravel(), Y.ravel())).reshape(X.shape)
h = 0.004
g = np.arange(-R, R + h, h)
X, Y = np.meshgrid(g, g)
d = dist(X, Y)
print("grid max", d.max(), " Lipschitz upper bound", d.max() + h / np.sqrt(2))
mf = ndimage.maximum_filter(d, size=7)
ii, jj = np.where((d == mf) & (d > 0.9))
out = []
for i, j in zip(ii, jj):
    cx, cy, hh = X[i, j], Y[i, j], h
    for it in range(30):
        gg = np.linspace(-2 * hh, 2 * hh, 41)
        XX, YY = np.meshgrid(cx + gg, cy + gg)
        dd = dist(XX, YY)
        k = np.unravel_index(dd.argmax(), dd.shape)
        cx, cy, hh = XX[k], YY[k], hh / 4
        if hh < 1e-13: break
    out.append((dd.max(), cx, cy))
out.sort(reverse=True)
seen = []
for v, cx, cy in out:
    if any(abs(cx - a) + abs(cy - b) < 1e-5 for _, a, b in seen): continue
    seen.append((v, cx, cy))
    # touching pieces
    idx = tree.query_ball_point([cx, cy], v + 1e-6)
    pcs = sorted(set(lab[idx]))
    desc = [f"{k}{S[k][0]}" + ("(tip)" if any(abs(allp[q]) < S[k][2] * (1 + 1e-7) for q in idx if lab[q] == k) else "") for k in pcs]
    if R - np.hypot(cx, cy) < v + 1e-6: desc.append("circle")
    print(f"r={v:.10f} |c|={np.hypot(cx,cy):.8f} arg={np.degrees(np.arctan2(cy,cx))%360:9.4f} touches {desc}")
