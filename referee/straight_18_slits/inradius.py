"""Independent inradius computation: global grid search over the whole disc
(no symmetry assumed) followed by local polishing of every local maximum."""
import numpy as np
from scipy.optimize import minimize
from geom import R, DEG, claimed_slits, dist_to_boundary

slits = claimed_slits()


def all_dists(p):
    """vector of distances from point p to each of the 19 boundary pieces"""
    out = [R - abs(p)]
    for a, t in slits:
        q = p * np.exp(-1j * a)
        x = min(max(q.real, t), R)
        out.append(np.hypot(q.real - x, q.imag))
    return np.array(out)


def polish(p0):
    """maximise min-distance by SLSQP on epigraph form: max r s.t. d_i(p) >= r"""
    x0 = np.array([p0.real, p0.imag, all_dists(p0).min()])
    cons = {'type': 'ineq', 'fun': lambda x: all_dists(x[0] + 1j * x[1]) - x[2]}
    res = minimize(lambda x: -x[2], x0, constraints=[cons], method='SLSQP',
                   options={'ftol': 1e-15, 'maxiter': 500})
    p = res.x[0] + 1j * res.x[1]
    return p, all_dists(p).min()


def main():
    h = 0.004
    xs = np.arange(-R, R + h, h)
    X, Y = np.meshgrid(xs, xs)
    P = X + 1j * Y
    D = np.empty(P.shape)
    for i in range(P.shape[0]):
        D[i], _ = dist_to_boundary(P[i], slits, R)
    D[np.abs(P) >= R] = -1
    print("grid", P.shape, "h =", h, " grid max d =", D.max(),
          "at", P.flat[D.argmax()])
    # Lipschitz bound: d is 1-Lipschitz so true sup <= gridmax + h/sqrt(2)
    print("Lipschitz upper bound on inradius from grid alone:", D.max() + h / np.sqrt(2))
    # local maxima on the grid (3x3) with d > 0.9
    c = D[1:-1, 1:-1]
    ismax = np.ones(c.shape, bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                ismax &= c >= D[1 + di:D.shape[0] - 1 + di, 1 + dj:D.shape[1] - 1 + dj]
    cand = P[1:-1, 1:-1][ismax & (c > 0.90)]
    print("number of grid local maxima with d>0.90:", cand.size)
    results = []
    for p0 in cand:
        p, r = polish(p0)
        results.append((r, p))
    # merge duplicates
    results.sort(key=lambda t: -t[0])
    uniq = []
    for r, p in results:
        if all(abs(p - q) > 1e-4 for _, q in uniq):
            uniq.append((r, p))
    names = ['circle'] + ['slit@%.5f deg tip %.5f' % (a / DEG, t) for a, t in slits]
    print("\nDistinct local maxima of the distance function (r, |p|, arg p deg, touching):")
    for r, p in uniq:
        d = all_dists(p)
        touch = [names[i] for i in np.where(d < r + 1e-7)[0]]
        print("  r=%.12f  |p|=%.9f arg=%9.5f  touches: %s" %
              (r, abs(p), np.angle(p) / DEG % 360, "; ".join(touch)))
    print("\nINRADIUS =", "%.12f" % uniq[0][0])
    # excess of the second family above/below 1
    print("values minus 1:", ["%.3e" % (r - 1) for r, _ in uniq])


if __name__ == "__main__":
    main()
