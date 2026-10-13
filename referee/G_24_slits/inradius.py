"""Step 2: inradius of D, whole disc, no symmetry.
 (a) branch-and-bound upper bound: d(p)=dist(p, boundary) is 1-Lipschitz; on a square of half-diagonal
     delta centred c, sup d <= d(c)+delta.  d(c) is computed with slits replaced by dense samples, which
     can only OVERestimate d -> the bound stays valid.
 (b) local maxima: grid h=0.004 -> KD-tree nested zoom -> Newton on equidistance system with exact
     nearest points on the exact curves (Brent polish) -> verify inscribed, list contacts.
"""
import sys, json
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import minimize_scalar, fsolve
from scipy import ndimage
from dom import R_G as R, slits_G, zslit

S = slits_G()
NS = 100001
XS = np.linspace(-1, 0, NS)
SAMP = [zslit(XS, R, *t[1:]) for t in S]
allp = np.concatenate(SAMP)
tree = cKDTree(np.c_[allp.real, allp.imag])
print("max sample spacing:", max(np.abs(np.diff(p)).max() for p in SAMP))


def dsamp(cx, cy):
    d, _ = tree.query(np.c_[cx, cy], workers=-1)
    return np.minimum(d, R - np.hypot(cx, cy))


def dist_slit(p, k):
    """exact distance to slit k: all sampled local minima polished by bounded Brent."""
    par = S[k][1:]
    d = np.abs(SAMP[k] - p)
    loc = list(np.where((d[1:-1] <= d[:-2]) & (d[1:-1] <= d[2:]))[0] + 1) + [0, NS - 1]
    dm = d.min()
    best = (np.inf, None)
    f = lambda x: abs(zslit(x, R, *par) - p)
    for i in loc:
        if d[i] > dm + 1e-3: continue
        lo, hi = XS[max(i - 2, 0)], XS[min(i + 2, NS - 1)]
        r = minimize_scalar(f, bounds=(lo, hi), method="bounded", options=dict(xatol=1e-16))
        for xx in (r.x, lo, hi):
            v = f(xx)
            if v < best[0]: best = (v, xx)
    return best


def pieces(p):
    return [(R - abs(p), None)] + [dist_slit(p, k) for k in range(len(S))]


def bnb(target, stop):
    side = 0.05
    g = np.arange(-R - side, R + side, side) + side / 2
    CX, CY = [a.ravel() for a in np.meshgrid(g, g)]
    ub = None
    while True:
        delta = side / np.sqrt(2)
        d = dsamp(CX, CY)
        keep = d + delta > target
        ub = (d + delta).max()
        CX, CY = CX[keep], CY[keep]
        print(f"  side={side:.3e} boxes={len(d):8d} survivors={keep.sum():8d} upper bound on R(D) so far={ub:.12f}", flush=True)
        if not keep.any() or side / 2 < stop: break
        q = side / 4
        CX = np.concatenate([CX - q, CX + q, CX - q, CX + q]); CY = np.concatenate([CY - q, CY - q, CY + q, CY + q])
        side /= 2
    return ub, CX + 1j * CY


def polish(p):
    ds = pieces(p); dmin = min(d for d, _ in ds)
    ks = [k for k, (dd, x) in enumerate(ds) if dd < dmin + 1e-5]
    q = p
    if len(ks) >= 3:
        best = None
        # try all triples of near-contacts, keep the solution that is inscribed with largest r
        from itertools import combinations
        for tri in combinations(ks, 3):
            def F(u):
                z = u[0] + 1j * u[1]
                v = [(R - abs(z)) if k == 0 else dist_slit(z, k - 1)[0] for k in tri]
                return [v[0] - u[2], v[1] - u[2], v[2] - u[2]]
            sol, info, ier, msg = fsolve(F, [p.real, p.imag, dmin], xtol=1e-15, full_output=True)
            z = sol[0] + 1j * sol[1]
            if abs(z - p) > 1e-3: continue
            rr = min(d for d, _ in pieces(z))
            if best is None or rr > best[0]: best = (rr, z)
        if best: q = best[1]
    ds = pieces(q); rr = min(d for d, _ in ds)
    desc, cps = [], []
    for k, (dd, x) in enumerate(ds):
        if dd < rr + 1e-8:
            if k == 0:
                desc.append("circle"); cps.append(q / abs(q) * R)
            else:
                desc.append(f"{k-1}{S[k-1][0]}:" + ("TIP" if x <= -1 + 1e-12 else f"x={x:.6f}"))
                cps.append(zslit(x, R, *S[k - 1][1:]))
    an = np.sort(np.angle(np.array(cps) - q)); gaps = np.diff(np.append(an, an[0] + 2 * np.pi))
    return rr, q, desc, np.degrees(gaps.max())


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "bnb":
        target, stop = float(sys.argv[2]), float(sys.argv[3])
        print(f"B&B: discarding boxes proven to have sup d <= {target}")
        ub, pts = bnb(target, stop)
        cl = []
        for p in pts:
            for c in cl:
                if abs(c[0] - p) < 2e-2: c.append(p); break
            else: cl.append([p])
        print(f"  final Lipschitz upper bound R(D) <= {ub:.12f}; {len(pts)} surviving boxes in {len(cl)} clusters")
        res = [polish(np.mean(c)) + (len(c), max(abs(np.array(c) - np.mean(c)))) for c in cl]
        res.sort(key=lambda t: -t[0])
        for rr, q, desc, gap, m, sp in res:
            print(f"   cluster({m} boxes, spread {sp:.1e}) exact max r={rr:.13f} (r-1={rr-1:+.2e}) |c|={abs(q):.9f} arg={np.degrees(np.angle(q))%360:9.5f} maxgap={gap:5.1f} {desc}")
        if res: print("  LOWER bound (exact inscribed disc): R(D) >= %.13f" % res[0][0])
        else: print("  no survivors: R(D) <= target proven (floating point, not interval arithmetic)")
    else:
        h = 0.004
        g = np.arange(-R, R + h, h)
        X, Y = np.meshgrid(g, g)
        d = dsamp(X.ravel(), Y.ravel()).reshape(X.shape)
        print("grid", X.shape, "grid max", d.max(), " crude Lipschitz upper bound", d.max() + h / np.sqrt(2))
        thr = float(sys.argv[2]) if len(sys.argv) > 2 else 0.97
        mf = ndimage.maximum_filter(d, size=9)
        ii, jj = np.where((d == mf) & (d > thr))
        print("grid local maxima above", thr, ":", len(ii))
        starts = []
        for i, j in zip(ii, jj):
            cx, cy, hh = X[i, j], Y[i, j], h
            for it in range(30):
                gg = np.linspace(-2 * hh, 2 * hh, 41)
                XX, YY = np.meshgrid(cx + gg, cy + gg)
                dd = dsamp(XX.ravel(), YY.ravel()).reshape(XX.shape)
                k = np.unravel_index(dd.argmax(), dd.shape)
                cx, cy, hh = XX[k], YY[k], hh / 4
                if hh < 1e-11: break
            if all(abs(cx + 1j * cy - s) > 1e-4 for s in starts): starts.append(cx + 1j * cy)
        res = [polish(p) for p in starts]
        res.sort(key=lambda t: -t[0])
        out = []
        for rr, q, desc, gap in res:
            flag = "" if (gap < 180 + 1e-6 and len(desc) >= 3) else "  [not strict max: ridge/degenerate]"
            print(f"r={rr:.13f} (r-1={rr-1:+.3e}) |c|={abs(q):.10f} arg={np.degrees(np.angle(q))%360:9.5f} maxgap={gap:6.1f} {desc}{flag}")
            out.append(dict(r=rr, x=q.real, y=q.imag, touch=desc))
        json.dump(out, open("maxima.json", "w"), indent=1)
        print("largest exact inscribed radius: %.15f" % res[0][0])
