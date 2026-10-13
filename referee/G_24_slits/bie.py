"""Step 3: conformal radius of (disc radius R minus slits) at 0.  z-plane BIE, no symmetry used.

    G0(z,zeta) = -log|z-zeta| + log|R^2 - conj(zeta) z| - log R     (Dirichlet Green fn of disc)
    g(z) = G0(z,0) - int_slits G0(z,zeta) mu(zeta)|dzeta|,  g=0 on slits
    log r(D) = log R - int_slits log(R/|zeta|) mu |dzeta|.

Unknown: sigma(t) = mu |dzeta/dt| per slit, slit parameter x = -1 + t^power (tip at t=0) so that the
inverse-sqrt tip singularity is absorbed.  n panels in t (edges optionally graded), p Gauss-Legendre
nodes per panel, Nystrom collocation at the nodes.
 * same slit: product integration, panel split at the (clamped) target parameter, dyadic refinement
   (L levels) toward the target, density interpolated by the panel's Legendre interpolant.
 * other slits: every source panel is oversampled into `m` sub-panels of p nodes with the density
   interpolated (m=1 -> plain Gauss-Legendre).  m is a convergence knob for the near interactions
   between close slits.  Optionally `adapt` : targets closer than `adapt`*panel length to a source
   panel are treated with dyadic refinement toward the nearest parameter instead (near-singular safe).
Adapted from the idea of referee_bloch_curved/bie.py but rewritten (b2 term, oversampling, adaptivity).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss, legvander
from dom import zslit


def G0(z, zeta, R):
    return -np.log(np.abs(z - zeta)) + np.log(np.abs(R * R - np.conj(zeta) * z)) - np.log(R)


def make_edges(n, grade_tip=1.0, grade_circ=1.0):
    """edges on [0,1]; grade>1 clusters panels at that end (algebraically)."""
    u = np.linspace(0, 1, n + 1)
    if grade_tip != 1.0 or grade_circ != 1.0:
        # map via regularised incomplete-beta like blend
        e = u ** grade_tip
        e = 1 - (1 - e) ** grade_circ
        return e
    return u


def conformal_radius(R, slits, n=12, p=16, L=40, power=2, m=1, edges=None, adapt=0.0, return_all=False):
    """slits: list of (alpha, tip, b1, b2)."""
    gx, gw = leggauss(p)
    if edges is None:
        edges = np.linspace(0, 1, n + 1)
    n = len(edges) - 1
    mid = 0.5 * (edges[1:] + edges[:-1]); hl = 0.5 * (edges[1:] - edges[:-1])
    t = (mid[:, None] + hl[:, None] * gx[None, :]).ravel()
    w = (hl[:, None] * gw[None, :]).ravel()
    M = n * p
    ns = len(slits)
    N = ns * M
    xpar = lambda tt: -1 + tt ** power
    Z = np.concatenate([zslit(xpar(t), R, *s) for s in slits])
    A = np.empty((N, N))
    Cinv = np.linalg.inv(legvander(gx, p - 1))          # nodal values -> Legendre coeffs
    # oversampled rule for cross interactions
    sub = np.linspace(-1, 1, m + 1)
    xi_f = ((sub[1:] + sub[:-1])[:, None] / 2 + (sub[1:] - sub[:-1])[:, None] / 2 * gx[None, :]).ravel()   # in [-1,1]
    wi_f = ((sub[1:] - sub[:-1])[:, None] / 2 * gw[None, :]).ravel()
    Pint = legvander(xi_f, p - 1) @ Cinv                 # (m p) x p interpolation
    t_f = (mid[:, None] + hl[:, None] * xi_f[None, :]).ravel()
    w_f = (hl[:, None] * wi_f[None, :]).ravel()
    # dyadic rule on [0,1] refined toward 1
    lo = 1 - 0.5 ** np.arange(L + 1)
    hi = np.append(1 - 0.5 ** np.arange(1, L + 1), 1.0)
    u = ((lo + hi)[:, None] / 2 + (hi - lo)[:, None] / 2 * gx[None, :]).ravel()
    uw = ((hi - lo)[:, None] / 2 * gw[None, :]).ravel()

    def product_block(zi, s, P, ts):
        """integrate G0(zi, .) * Lagrange basis of panel P of slit s; split at ts (clamped)."""
        a, b = edges[P], edges[P + 1]
        ts = np.clip(ts, a, b)
        sL = a + (ts - a)[:, None] * u[None, :]; wL = (ts - a)[:, None] * uw[None, :]
        sR = b - (b - ts)[:, None] * u[None, :]; wR = (b - ts)[:, None] * uw[None, :]
        sq = np.concatenate([sL, sR], axis=1); wq = np.concatenate([wL, wR], axis=1)
        zeta = zslit(xpar(sq), R, *s)
        with np.errstate(divide="ignore", invalid="ignore"):
            K = G0(zi[:, None], zeta, R)
        K = np.where(wq > 0, K, 0.0)
        K[~np.isfinite(K)] = 0.0
        V = legvander((2 * sq - a - b) / (b - a), p - 1)
        return np.einsum("mq,mqj->mj", wq * K, V) @ Cinv

    nadapt = 0
    for k, s in enumerate(slits):
        sl = slice(k * M, (k + 1) * M)
        # cross interactions: all targets, source slit k
        Zf = zslit(xpar(t_f), R, *s)
        with np.errstate(divide="ignore"):
            Kf = G0(Z[:, None], Zf[None, :], R) * w_f[None, :]
        A[:, sl] = (Kf.reshape(N, n, m * p) @ Pint).reshape(N, M)
        if adapt > 0:
            dmin = np.abs(Z[:, None] - Zf[None, :]).reshape(N, n, m * p)
            imin = dmin.argmin(axis=2); dmin = dmin.min(axis=2)
            ze = zslit(xpar(edges), R, *s); plen = np.abs(np.diff(ze))
            near = dmin < adapt * plen[None, :]
            near[sl, :] = False
            tfm = t_f.reshape(n, m * p)
            for P in range(n):
                idx = np.where(near[:, P])[0]
                if len(idx) == 0: continue
                nadapt += len(idx)
                # nearest parameter (refine by a few golden steps on dense sample)
                tt = np.linspace(edges[P], edges[P + 1], 2001)
                zz = zslit(xpar(tt), R, *s)
                ts = tt[np.abs(Z[idx][:, None] - zz[None, :]).argmin(axis=1)]
                A[idx, k * M + P * p: k * M + (P + 1) * p] = product_block(Z[idx], s, P, ts)
        # self interactions
        zi = Z[sl]
        for P in range(n):
            A[sl, k * M + P * p: k * M + (P + 1) * p] = product_block(zi, s, P, t)
    rhs = np.log(R / np.abs(Z))
    sig = np.linalg.solve(A, rhs)
    W = np.tile(w, ns)
    logr = np.log(R) - np.sum(sig * W * rhs)
    if return_all:
        return np.exp(logr), dict(t=t, w=w, Z=Z, sig=sig, edges=edges, A=A, nadapt=nadapt, power=power)
    return np.exp(logr)
