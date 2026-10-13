"""Step 3: conformal radius of (disc of radius R minus slits) at 0.

Method (z-plane, NOT the claimant's): single-layer potential on the slits with the Dirichlet
Green's function of the disc,
    G0(z,zeta) = -log|z-zeta| + log|R^2 - conj(zeta) z| - log R,
    g(z) = G0(z,0) - int_slits G0(z,zeta) mu(zeta) |dzeta|,   g = 0 on slits,
    log r(D) = log R - int_slits log(R/|zeta|) mu |dzeta|.
Each slit is parametrised by x = -1 + t^2, t in [0,1] (so that mu |dzeta/dt| is smooth at the
tip), cut into n panels in t, with p Gauss-Legendre nodes per panel (piecewise polynomial
density, collocation at the nodes).  Interactions between different slits use the plain
Gauss-Legendre rule; interactions of a target with every panel of its own slit use product
integration: the panel is split at the (clamped) target parameter and integrated with a
dyadically refined composite Gauss rule against the Lagrange basis.
No symmetry of the domain is used.
"""
import sys, time
import numpy as np
from numpy.polynomial.legendre import leggauss, legvander
from dom import zslit


def G0(z, zeta, R):
    return -np.log(np.abs(z - zeta)) + np.log(np.abs(R * R - np.conj(zeta) * z)) - np.log(R)


def conformal_radius(R, slits, n=12, p=16, L=44, power=2, verbose=False, edges=None):
    """slits: list of (alpha, tip, bend)."""
    gx, gw = leggauss(p)
    if edges is None:
        edges = np.linspace(0, 1, n + 1)
    n = len(edges) - 1
    mid = 0.5 * (edges[1:] + edges[:-1])
    hl = 0.5 * (edges[1:] - edges[:-1])
    t = (mid[:, None] + hl[:, None] * gx[None, :]).ravel()
    w = (hl[:, None] * gw[None, :]).ravel()
    M = n * p
    ns = len(slits)
    N = ns * M
    xpar = lambda tt: -1 + tt ** power
    Z = np.concatenate([zslit(xpar(t), R, *s) for s in slits])
    W = np.tile(w, ns)
    with np.errstate(divide="ignore", invalid="ignore"):
        A = G0(Z[:, None], Z[None, :], R) * W[None, :]
    Cinv = np.linalg.inv(legvander(gx, p - 1))
    # dyadic composite rule on [0,1] refined toward 1:  pieces [1-2^-l, 1-2^-(l+1)], l<L, + last
    lo = 1 - 0.5 ** np.arange(L + 1)
    hi = np.append(1 - 0.5 ** np.arange(1, L + 1), 1.0)
    u = ((lo + hi)[:, None] / 2 + (hi - lo)[:, None] / 2 * gx[None, :]).ravel()   # in [0,1]
    uw = ((hi - lo)[:, None] / 2 * gw[None, :]).ravel()
    for k, s in enumerate(slits):
        sl = slice(k * M, (k + 1) * M)
        zi = Z[sl]
        for P in range(n):
            a, b = edges[P], edges[P + 1]
            ts = np.clip(t, a, b)
            # left side [a, ts]: s = a + (ts-a) u ; right side [ts, b]: s = b - (b-ts) u
            sL = a + (ts - a)[:, None] * u[None, :]
            wL = (ts - a)[:, None] * uw[None, :]
            sR = b - (b - ts)[:, None] * u[None, :]
            wR = (b - ts)[:, None] * uw[None, :]
            sq = np.concatenate([sL, sR], axis=1)
            wq = np.concatenate([wL, wR], axis=1)
            zeta = zslit(xpar(sq), R, *s)
            with np.errstate(divide="ignore", invalid="ignore"):
                K = G0(zi[:, None], zeta, R)
            K = np.where(wq > 0, K, 0.0)
            K[~np.isfinite(K)] = 0.0
            xi = (2 * sq - a - b) / (b - a)
            V = legvander(xi, p - 1)
            blk = np.einsum("mq,mqj->mj", wq * K, V) @ Cinv
            A[sl, k * M + P * p: k * M + (P + 1) * p] = blk
    rhs = np.log(R / np.abs(Z))
    mu = np.linalg.solve(A, rhs)
    logr = np.log(R) - np.sum(mu * W * np.log(R / np.abs(Z)))
    if verbose:
        print("   N=%d cond~%.2e total charge %.12f" % (N, np.linalg.cond(A) if N < 3000 else -1, np.sum(mu * W)))
    return np.exp(logr), (t, w, Z, mu)


if __name__ == "__main__":
    # calibration 1: unit disc minus [s,1]
    for s in (0.5, 0.1, 0.9):
        ex = 4 * s / (1 + s) ** 2
        for n in (2, 4, 8):
            r, _ = conformal_radius(1.0, [(0.0, s, 0.0)], n=n)
            print(f"cal1 s={s} n={n}: r={r:.15f} exact={ex:.15f} err={r-ex:+.2e}")
