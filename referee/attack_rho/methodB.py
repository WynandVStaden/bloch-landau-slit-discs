"""Method B: panel Nystrom BIE with the DISC Green function on the physical slits only.

 * unknown: charge density on each slit in the tip-adapted variable tau (s = -a + a tau^2, removes the
   inverse square-root tip singularity);  dmu = g(tau) dtau,  tau in [0,1]
 * kernel:  G1(x,y) = log|1 - x conj(y)| - log|x - y|   (unit disc Green function; no reflected arcs,
   no unknown constant)
 * discretisation: composite Gauss-Legendre panels, product integration for the three (near-)singular
   logarithms  log|tau'-tau|, log(tau'+tau), log(tau*-tau')  (tau* = sqrt(2-tau^2) = mirror image)
 * variant B0: z-plane directly, all 5 classes with D3 symmetry
 * variant B1: zeta=(z/R)^3, then the two radial slit classes are removed EXACTLY by the explicit map
       F = Jinv o affine o J  of  (unit disc minus [p,1) minus (-1,-q])  onto the unit disc,
   leaving 6 curved slits (3 classes x conj).   rho(D) = R (rho_eta * 2/(A+B))^(1/3).
"""
import sys
import numpy as np
import mpmath as mp
from geom import full_problem, R0, CLASSES0

# ---------------------------------------------------------------- quadrature tools


def gauss_log_rule(n):
    """Gauss rule for int_0^1 (-log v) f(v) dv (Golub-Welsch from exact moments, in high precision)."""
    old = mp.mp.dps
    mp.mp.dps = 120
    H = mp.matrix(n + 1, n + 1)
    for i in range(n + 1):
        for j in range(n + 1):
            H[i, j] = mp.mpf(1) / (i + j + 1)**2
    Rm = mp.cholesky(H).T          # upper triangular, H = R^T R
    al = [Rm[0, 1] / Rm[0, 0]] + [Rm[j, j + 1] / Rm[j, j] - Rm[j - 1, j] / Rm[j - 1, j - 1] for j in range(1, n)]
    be = [Rm[j + 1, j + 1] / Rm[j, j] for j in range(n - 1)]
    J = mp.matrix(n, n)
    for j in range(n):
        J[j, j] = al[j]
    for j in range(n - 1):
        J[j, j + 1] = be[j]
        J[j + 1, j] = be[j]
    E, Q = mp.eigsy(J)
    x = np.array([float(E[i]) for i in range(n)])
    w = np.array([float(Q[0, i]**2) for i in range(n)])     # m0 = 1
    mp.mp.dps = old
    idx = np.argsort(x)
    return x[idx], w[idx]


_cache = {}


def gl01(n):
    if n not in _cache:
        x, w = np.polynomial.legendre.leggauss(n)
        _cache[n] = ((x + 1) / 2, w / 2)
    return _cache[n]


def glog(n):
    k = ('log', n)
    if k not in _cache:
        _cache[k] = gauss_log_rule(n)
    return _cache[k]


NQ = 20   # order of auxiliary rules


def log_rule(x, c, d):
    """nodes u, weights om with int_c^d log|x-y| p(y) dy ~ sum om p(u) for smooth p (x anywhere on the line)."""
    gx, gw = gl01(NQ)
    us, oms = [], []
    if c < x < d:
        lx, lw = glog(NQ)
        for L, sg in ((x - c, -1.0), (d - x, 1.0)):
            us += [x + sg * L * gx, x + sg * L * lx]
            oms += [L * np.log(L) * gw, -L * lw]
    else:
        if x <= c:
            e, far, sg = c, d, 1.0
        else:
            e, far, sg = d, c, -1.0
        delta = abs(x - e)
        tot = abs(far - e)
        lo = 0.0
        ln = max(delta, 1e-300)
        while lo < tot:
            hi = min(lo + ln, tot)
            y = e + sg * (lo + (hi - lo) * gx)
            us.append(y)
            oms.append((hi - lo) * gw * np.log(np.abs(x - y)))
            lo = hi
            ln *= 2
    return np.concatenate(us), np.concatenate(oms)


def bary_matrix(xn, xe):
    """Lagrange interpolation matrix from nodes xn to points xe."""
    w = np.array([1.0 / np.prod(xn[k] - np.delete(xn, k)) for k in range(len(xn))])
    D = xe[:, None] - xn[None, :]
    D[D == 0] = 1e-300
    M = w[None, :] / D
    return M / M.sum(axis=1, keepdims=True)


# ---------------------------------------------------------------- arcs

class Arc:
    """slit class in the eta-plane: eta(tau), d eta/d sigma (sigma = tau^2), images = list of callables."""

    def __init__(self, a, theta, b1, b2, images, premap=None, power=1):
        self.a, self.theta, self.b1, self.b2 = a, theta, b1, b2
        self.images = images
        self.premap = premap
        self.power = power

    def eta(self, tau, deriv=False):
        s = -self.a + self.a * tau**2
        v = self.theta + self.b1 * s**2 + self.b2 * s**4
        z = np.exp(self.power * (s + 1j * v))
        dz = z * self.power * (1 + 1j * (2 * self.b1 * s + 4 * self.b2 * s**3)) * self.a   # d/d sigma
        if self.premap is None:
            return (z, dz) if deriv else z
        e, de = self.premap(z)
        return (e, de * dz) if deriv else e


def two_slit_map(p, q):
    """F: unit disc minus [p,1) minus (-1,-q]  ->  unit disc, F(0)=0, F'(0)=(A+B)/2.  returns (F, F')"""
    A, B = (p + 1 / p) / 2, (q + 1 / q) / 2

    def F(zeta):
        Jz = (zeta + 1 / zeta) / 2
        u = (2 * Jz - (A - B)) / (A + B)
        r = np.sqrt(u * u - 1)
        e1, e2 = u - r, u + r
        e = np.where(np.abs(e1) < np.abs(e2), e1, e2)
        r = (1 / e - e) / 2
        de = (-e / r) * (1 - 1 / zeta**2) / (A + B)
        return e, de
    return F, (A + B) / 2


# ---------------------------------------------------------------- solver

def panels(P, grade_base=0, grade_tip=0):
    """panel break points on [0,1]; optional dyadic refinement toward tau=1 (base) / tau=0 (tip)."""
    br = list(np.linspace(0, 1, P + 1))
    h = 1.0 / P
    for k in range(1, grade_base + 1):
        br.append(1 - h / 2**k)
    for k in range(1, grade_tip + 1):
        br.append(h / 2**k)
    return np.array(sorted(br))


def solve(arcs, P, q, grade_base=0, grade_tip=0):
    gx, gw = gl01(q)
    br = panels(P, grade_base, grade_tip)
    npan = len(br) - 1
    tau = np.concatenate([br[k] + (br[k + 1] - br[k]) * gx for k in range(npan)])
    wt = np.concatenate([(br[k + 1] - br[k]) * gw for k in range(npan)])
    n1 = len(tau)
    nc = len(arcs)
    A = np.zeros((nc * n1, nc * n1))
    rhs = np.zeros(nc * n1)
    E = [arc.eta(tau, deriv=True) for arc in arcs]
    # singular weights (geometry independent): Wm[i,j] for -log|t'-t| - log(t'+t) + log(t*-t')
    Ws = np.zeros((n1, n1))
    tstar = np.sqrt(2 - tau**2)
    for k in range(npan):
        c, d = br[k], br[k + 1]
        xn = tau[k * q:(k + 1) * q]
        for i in range(n1):
            row = np.zeros(q)
            for xs, sgn in ((tau[i], -1.0), (-tau[i], -1.0), (tstar[i], 1.0)):
                u, om = log_rule(xs, c, d)
                row += sgn * (om @ bary_matrix(xn, u))
            Ws[i, k * q:(k + 1) * q] = row
    for i, ai in enumerate(arcs):
        ei, dei = E[i]
        for j, aj in enumerate(arcs):
            ej, dej = E[j]
            blk = np.zeros((n1, n1))
            for im, img in enumerate(aj.images):
                y = img(ej)
                X, Y = ei[:, None], y[None, :]
                if i == j and im == 0:
                    sg_i, sg_j = tau[:, None]**2, tau[None, :]**2
                    dsg = sg_i - sg_j
                    np.fill_diagonal(dsg, 1.0)
                    num = (X - Y) / dsg
                    num[np.diag_indices(n1)] = dei
                    S = -np.log(np.abs(num)) + np.log(np.abs(1 - X * np.conj(Y)) / (tstar[:, None] - tau[None, :]))
                    blk += Ws + S * wt[None, :]
                else:
                    blk += (np.log(np.abs(1 - X * np.conj(Y))) - np.log(np.abs(X - Y))) * wt[None, :]
            A[i * n1:(i + 1) * n1, j * n1:(j + 1) * n1] = blk
        rhs[i * n1:(i + 1) * n1] = -np.log(np.abs(ei))
    g = np.linalg.solve(A, rhs).reshape(nc, n1)
    logrho = 0.0
    hm = []
    for i, ai in enumerate(arcs):
        logrho -= len(ai.images) * np.sum(wt * g[i] * (-np.log(np.abs(E[i][0]))))
        hm.append(len(ai.images) * np.sum(wt * g[i]))
    return np.exp(logrho), np.array(hm), dict(tau=tau, wt=wt, g=g, cond=np.linalg.cond(A), E=E)


def rot(th):
    return lambda z: z * np.exp(1j * th)


def rotconj(th):
    return lambda z: np.conj(z) * np.exp(1j * th)


def arcs_B0(R=R0, classes=CLASSES0, nrot=3):
    arcs = []
    for i, (th, t, b1, b2) in enumerate(classes):
        imgs = [rot(2 * np.pi * r / nrot) for r in range(nrot)]
        if not (b1 == 0 and b2 == 0 and (abs(np.sin(nrot * th)) < 1e-12)):
            imgs += [rotconj(2 * np.pi * r / nrot) for r in range(nrot)]
        arcs.append(Arc(np.log(R / t), th, b1, b2, imgs))
    return arcs, (lambda rho_eta: R * rho_eta)


def arcs_B1(R=R0, classes=CLASSES0):
    p = (classes[0][1] / R)**3
    qq = (classes[1][1] / R)**3
    F, F0 = two_slit_map(p, qq)
    arcs = []
    for (th, t, b1, b2) in classes[2:]:
        arcs.append(Arc(np.log(R / t), th, b1, b2, [lambda z: z, np.conj], premap=F, power=3))
    return arcs, (lambda rho_eta: R * (rho_eta / F0)**(1.0 / 3))


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'B1'
    arcs, back = arcs_B0() if which == 'B0' else arcs_B1()
    prev = None
    for (P, q, gb, gt) in [(2, 8, 0, 0), (4, 8, 0, 0), (4, 12, 0, 0), (8, 12, 0, 0), (8, 16, 0, 0), (12, 16, 0, 0),
                           (8, 16, 3, 0), (12, 16, 4, 2), (16, 16, 6, 0)]:
        rho_eta, hm, info = solve(arcs, P, q, gb, gt)
        rho = back(rho_eta)
        d = '' if prev is None else f'{rho - prev:+.3e}'
        print(f'{which} P={P:2d} q={q:2d} gb={gb} gt={gt} n/arc={len(info["tau"]):4d} rho={rho:.15f} diff={d:>11s} '
              f'cond={info["cond"]:.1e}', flush=True)
        prev = rho
    print('harmonic measure per class:', hm, ' sum:', hm.sum())
