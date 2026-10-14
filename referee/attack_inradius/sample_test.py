# Referee: random-point attack on the stored leaves. Every sampled point of (sector ∩ disc) must lie in a leaf of
# kind (c) whose witness q (a genuine boundary point, recomputed with the referee's slit table) satisfies |p-q| <= 1+1e-9.
import numpy as np, pickle, time
from geom import *
z = np.load("stored_leaves_G.npz"); C, PARS = z["cells"], z["pars"]
LAM = 1 + 1e-9
levels = np.unique(C[:, 0]); keys = {}
for lv in levels:
    m = np.nonzero(C[:, 0] == lv)[0]; k = (C[m, 1].astype(np.uint64) << np.uint64(32)) | C[m, 2].astype(np.uint64)
    o = np.argsort(k); keys[int(lv)] = (k[o], m[o])
A = np.array([s['a'] for s in SL_F]); TH = np.array([s['th'] for s in SL_F]); B1 = np.array([s['b1'] for s in SL_F]); B2 = np.array([s['b2'] for s in SL_F])
def witness(rows):
    kind, j, par = C[rows, 4], C[rows, 5], PARS[rows]
    s = np.where(kind == 1, -A[j], par)
    assert np.all((kind == 2) | ((s <= 0) & (s >= -A[j])))
    q = R_F*np.exp(s + 1j*(TH[j] + B1[j]*s**2 + B2[j]*s**4))
    return np.where(kind == 2, R_F*np.exp(1j*par), q)
def locate(p):
    x, y = p.real, p.imag; leaf = np.full(len(p), -1); nfound = np.zeros(len(p), int)
    for lv, (ks, rows) in keys.items():
        ix = np.floor(x*2.0**(lv - 3)).astype(np.int64); iy = np.floor(y*2.0**(lv - 3)).astype(np.int64)
        ok = (ix >= 0) & (iy >= 0) & (ix < 2**lv) & (iy < 2**lv)
        k = (ix.astype(np.uint64) << np.uint64(32)) | iy.astype(np.uint64)
        pos = np.searchsorted(ks, k); pos[pos >= len(ks)] = len(ks) - 1
        hit = ok & (ks[pos] == k)
        leaf[hit] = rows[pos[hit]]; nfound += hit
    return leaf, nfound
from quadcover_lib import dist
def fold(p):
    a = np.mod(np.angle(p), 2*np.pi/3); a = np.where(a > np.pi/3, 2*np.pi/3 - a, a); return np.abs(p)*np.exp(1j*a)
def test(p, label):
    p = p[(np.abs(p) < R_F)]
    ang = np.angle(p); p = p[(ang >= 0) & (ang <= np.pi/3)]
    leaf, nf = locate(p)
    code = C[leaf, 3]
    bad_cover = (nf != 1); bad_code = (code != 2) & ~bad_cover
    q = witness(leaf); dq = np.abs(p - q); bad_w = (dq > LAM) & (code == 2)
    d = dist(p)
    print("%-44s n=%8d  not in exactly one leaf: %d  in (a)/(b) leaf: %d  witness farther than 1+1e-9: %d   max|p-q|-1=%+.3e  max d(p)-1=%+.3e  min(|p-q|-d(p))=%+.2e"
          % (label, len(p), bad_cover.sum(), bad_code.sum(), bad_w.sum(), dq[code == 2].max() - 1, d.max() - 1, (dq - d)[code == 2].min()), flush=True)
    if bad_code.sum():
        for pp, lf in list(zip(p[bad_code], leaf[bad_code]))[:5]: print("    point", pp, "abs", abs(pp), "arg", np.angle(pp), "leaf", C[lf])
    return bad_cover.sum() + bad_code.sum() + bad_w.sum()
rng = np.random.default_rng(12345); tot = 0; t0 = time.time()
for rep in range(4):
    n = 1_500_000; p = R_F*np.sqrt(rng.random(n))*np.exp(1j*np.pi/3*rng.random(n)); tot += test(p, "uniform in sector, batch %d" % rep)
# near the sector edges and the circle
n = 500_000
tot += test(R_F*np.sqrt(rng.random(n))*np.exp(1j*1e-9*rng.random(n)), "edge arg ~ 0+")
tot += test(R_F*np.sqrt(rng.random(n))*np.exp(1j*(np.pi/3 - 1e-9*rng.random(n))), "edge arg ~ pi/3-")
tot += test(R_F*(1 - 1e-9*rng.random(n))*np.exp(1j*np.pi/3*rng.random(n)), "near the circle")
tot += test(10.0**(-12*rng.random(n))*np.exp(1j*np.pi/3*rng.random(n)), "near the origin (log-uniform radius)")
mx = pickle.load(open("maxima.pkl", "rb"))
cs = np.array([complex(o['c']) for o in mx])
for k in range(1, 14):
    sc = 10.0**(-k); n = 60000
    pts = (cs[:, None] + sc*(rng.standard_normal((len(cs), n)) + 1j*rng.standard_normal((len(cs), n)))).ravel()
    tot += test(fold(pts), "all 22 maxima (folded), gaussian scale 1e-%d" % k)
# slit neighbourhoods (channels): points at distance <= 0.08 from slit points in the sector
for k in range(24):
    sl = SL_F[k]; s = -sl['a']*rng.random(40000); g = gam_f(sl, s) + 0.08*rng.random(40000)*np.exp(2j*np.pi*rng.random(40000))
    tot += test(fold(g), "around slit %d (folded)" % k) if k in (0, 3, 6, 12, 18) else 0
print("TOTAL FAILURES:", tot, " time %.0f s" % (time.time() - t0))
