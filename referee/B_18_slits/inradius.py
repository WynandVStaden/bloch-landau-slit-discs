"""Step 2: independent inradius, no symmetry assumed.
Stage 1: global grid (step h) with KD-tree distance to densely sampled slits + exact circle.
Stage 2: every grid local maximum above a threshold is refined with the exact distance
         function (nearest point on each exact curve found by sample + bounded Brent polish).
Stage 3: for the refined maxima, identify touching pieces and re-solve the equidistance
         system by Newton (fsolve) as a cross-check.
"""
import json, sys
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import minimize_scalar, fsolve
from scipy import ndimage
from dom import R_CLAIM as R, claimed_slits, zslit

S = claimed_slits()
NS = 4001
XS = np.linspace(-1, 0, NS)
SAMP = [zslit(XS, R, al, tip, b) for (c, al, tip, b) in S]


def dist_slit(p, k):
    """exact distance from complex p to slit k; returns (d, x*)"""
    c, al, tip, b = S[k]
    d = np.abs(SAMP[k] - p)
    # candidate local minima of the sampled distance
    idx = [d.argmin()]
    loc = np.where((d[1:-1] <= d[:-2]) & (d[1:-1] <= d[2:]))[0] + 1
    for i in loc:
        if d[i] < d[idx[0]] + 1e-3 and i != idx[0]:
            idx.append(i)
    best = (np.inf, None)
    for i in idx:
        lo = XS[max(i - 1, 0)]
        hi = XS[min(i + 1, NS - 1)]
        f = lambda x: abs(zslit(x, R, al, tip, b) - p)
        r = minimize_scalar(f, bounds=(lo, hi), method="bounded", options=dict(xatol=1e-15))
        for xx in (r.x, lo, hi):
            v = f(xx)
            if v < best[0]:
                best = (v, xx)
    return best


def dist_all(p):
    ds = [(R - abs(p), None)] + [dist_slit(p, k) for k in range(len(S))]
    return ds


def dexact(p):
    return min(d for d, _ in dist_all(p))


def refine(p0, h0):
    p = p0
    best = dexact(p)
    h = h0
    while h > 1e-14:
        improved = False
        g = np.arange(-3, 4) * h
        for dx in g:
            for dy in g:
                q = p0 + dx + 1j * dy
                v = dexact(q)
                if v > best:
                    best, p, improved = v, q, True
        p0 = p
        if not improved:
            h *= 0.5
    return p, best


def main():
    h = 0.004
    allp = np.concatenate(SAMP)
    tree = cKDTree(np.c_[allp.real, allp.imag])
    g = np.arange(-R, R + h, h)
    X, Y = np.meshgrid(g, g)
    d, _ = tree.query(np.c_[X.ravel(), Y.ravel()], workers=-1)
    d = d.reshape(X.shape)
    d = np.minimum(d, R - np.hypot(X, Y))
    d[np.hypot(X, Y) >= R] = -1
    print("grid", X.shape, "max grid dist", d.max())
    # Lipschitz: true sup <= grid max + h/sqrt2 (+ sampling error ~1e-7)
    print("crude rigorous upper bound from grid: ", d.max() + h / np.sqrt(2))
    mf = ndimage.maximum_filter(d, size=5)
    thr = 0.90
    ii, jj = np.where((d == mf) & (d > thr))
    print("grid local maxima above", thr, ":", len(ii))
    # every grid point above (best - h/sqrt2) must sit near some refined maximum -> check later
    res = []
    for i, j in zip(ii, jj):
        p, v = refine(X[i, j] + 1j * Y[i, j], h)
        res.append((v, p))
    res.sort(key=lambda t: -t[0])
    uniq = []
    for v, p in res:
        if all(abs(p - q) > 1e-6 for _, q in uniq):
            uniq.append((v, p))
    out = []
    print("\nrefined local maxima (r, centre mod, centre arg deg, touching pieces [idx class x*]):")
    for v, p in uniq:
        ds = dist_all(p)
        touch = [(k - 1, x) for k, (dd, x) in enumerate(ds) if dd < v + 1e-7]
        desc = []
        for k, x in touch:
            if k < 0:
                desc.append("circle")
            else:
                tag = "TIP" if x <= -1 + 1e-9 else ("x=%.6f" % x)
                desc.append(f"{k}{S[k][0]}:{tag}")
        print(f"  r={v:.12f}  |c|={abs(p):.9f} arg={np.degrees(np.angle(p)) % 360:9.4f}  {desc}")
        out.append(dict(r=v, x=p.real, y=p.imag, touch=desc))
    json.dump(out, open("inradius_maxima.json", "w"), indent=1)

    best = uniq[0][0]
    # coverage check: grid points with value > best - h (>= h/sqrt2 margin) all lie within
    # 0.3 of a refined maximum whose value we know exactly
    ci, cj = np.where(d > best - h)
    P = X[ci, cj] + 1j * Y[ci, cj]
    cen = np.array([p for _, p in uniq])
    dm = np.abs(P[:, None] - cen[None, :]).min(axis=1)
    print(f"\n{len(P)} grid points have d > best-h; max distance to a refined maximum: {dm.max():.4f}")
    print("best refined inradius: %.13f" % best)
    print("R(D) - 1 = %.3e" % (best - 1))

    # Newton cross-check on top maxima with exactly-3 smooth/tip contacts
    print("\nNewton equidistance cross-check for maxima with r > 0.999:")
    for v, p in uniq:
        if v < 0.999:
            break
        ds = dist_all(p)
        ks = [k for k, (dd, x) in enumerate(ds) if dd < v + 1e-7]

        def F(u):
            q = u[0] + 1j * u[1]
            vals = []
            for k in ks[:3]:
                vals.append(R - abs(q) if k == 0 else dist_slit(q, k - 1)[0])
            return [vals[0] - u[2], vals[1] - u[2], vals[2] - u[2]]
        if len(ks) >= 3:
            sol = fsolve(F, [p.real, p.imag, v], xtol=1e-15)
            print(f"  zoom r={v:.13f}  newton r={sol[2]:.13f}  pieces={[k-1 for k in ks]}")
        else:
            print(f"  zoom r={v:.13f}  only {len(ks)} contacts: {[k-1 for k in ks]}")


if __name__ == "__main__":
    main()
