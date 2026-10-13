"""Referee's independent conformal-radius computation.

Method: Schwarz-Christoffel-type parameter problem for radial-slit discs
(NOT an integral equation).

Let g: unit disc -> D = {|z|<R} minus radial slits, g(0)=0, g'(0)=c>0.
P(zeta) = zeta g'(zeta)/g(zeta) is analytic, zero-free in the disc, P(0)=1,
real on boundary arcs mapped into |z|=R and purely imaginary on arcs mapped
onto slit sides.  Hence P^2 is rational, and on zeta = e^{i phi}

    |P(phi)| = prod_j |sin((phi-alpha_j)/2)|
               / sqrt( prod_j |sin((phi-beta_j)/2) sin((phi-gamma_j)/2)| )

where alpha_j = prevertex of tip j, beta_j < alpha_j < gamma_j prevertices of
the two sides of the junction of slit j with the circle.

Side-length conditions:  int_{beta_j}^{alpha_j} |P| = int_{alpha_j}^{gamma_j} |P| = log(R/t_j)
Angle conditions:        int_{gamma_j}^{beta_{j+1}} |P| = theta_{j+1}-theta_j
Conformal radius (mean value property of log|g(zeta)/zeta|):
    log c = log R - (1/2pi) sum_j int_{beta_j}^{gamma_j} |phi-alpha_j| |P(phi)| dphi

All integrals: Gauss-Jacobi (exact treatment of the inverse-sqrt endpoint
singularities), on composite panels where needed.
"""
import numpy as np
from scipy.special import roots_jacobi, roots_legendre
from scipy.optimize import least_squares

TWO_PI = 2 * np.pi


class SlitDisc:
    def __init__(self, R, thetas, tips, n=64, npan=1):
        """thetas (rad, increasing in [0,2pi)), tips: slit data. n: nodes/panel."""
        self.R = R
        self.th = np.asarray(thetas, float)
        self.t = np.asarray(tips, float)
        self.N = len(self.th)
        self.n = n
        self.npan = npan
        # one-sided inverse sqrt weight at x=+1 : (1-x)^(-1/2)
        self.xj, self.wj = roots_jacobi(n, -0.5, 0.0)
        # both ends: Gauss-Chebyshev
        self.xc, self.wc = roots_jacobi(n, -0.5, -0.5)
        self.xl, self.wl = roots_legendre(n)

    # ---- |P| with selected singular factors removed analytically ----------
    def logabsP(self, phi, al, be, ga, skip=()):
        """log|P(phi)|, omitting the factors |sin((phi-s)/2)|^(-1/2) for the
        singular prevertices s listed in skip (those are put in the weight)."""
        phi = np.asarray(phi, float)[:, None]
        v = np.sum(np.log(np.abs(np.sin((phi - al[None, :]) / 2))), axis=1)
        sing = np.concatenate([be, ga])
        if len(skip):
            keep = np.ones(sing.size, bool)
            for s in skip:
                keep &= np.abs(np.sin((sing - s) / 2)) > 1e-14   # compare mod 2pi
            sing = sing[keep]
        v -= 0.5 * np.sum(np.log(np.abs(np.sin((phi - sing[None, :]) / 2))), axis=1)
        return v

    def int_side(self, a, b, al, be, ga, weightfun=None):
        """integral of |P| (times weightfun) over [a,b] (a may be > b) where b is
        a singular prevertex (inverse sqrt) and a is a tip prevertex (regular).
        Composite: (npan-1) Gauss-Legendre panels graded towards b + one
        Gauss-Jacobi panel adjacent to b."""
        L = b - a
        # geometric panel edges in fraction of the interval, graded towards b
        if self.npan == 1:
            edges = np.array([0.0, 1.0])
        else:
            edges = 1.0 - 0.5 ** np.arange(self.npan)
            edges = np.concatenate([edges, [1.0]])
        tot = 0.0
        for k in range(len(edges) - 1):
            p0 = a + L * edges[k]
            p1 = a + L * edges[k + 1]
            last = (k == len(edges) - 2)
            if last:
                x, w = self.xj, self.wj
                phi = 0.5 * (p0 + p1) + 0.5 * (p1 - p0) * x
                # singular factor |sin((phi-b)/2)|^(-1/2); weight is (1-x)^(-1/2)
                # with |phi-b| = |p1-p0|/2 * (1-x)
                f = np.exp(self.logabsP(phi, al, be, ga, skip=(b,)))
                d = np.abs(phi - b)
                corr = np.sqrt((d / 2) / np.abs(np.sin(d / 2)))  # smooth
                f = f * corr * (np.abs(p1 - p0) / 4) ** (-0.5)
            else:
                x, w = self.xl, self.wl
                phi = 0.5 * (p0 + p1) + 0.5 * (p1 - p0) * x
                f = np.exp(self.logabsP(phi, al, be, ga))
            if weightfun is not None:
                f = f * weightfun(phi)
            tot += np.sum(w * f) * abs(p1 - p0) / 2
        return tot

    def int_arc(self, a, b, al, be, ga):
        """integral of |P| over [a,b], both endpoints singular (inverse sqrt).
        Split at the midpoint into two one-sided integrals."""
        m = 0.5 * (a + b)
        return self.int_side(m, a, al, be, ga) + self.int_side(m, b, al, be, ga)

    # ---- residuals ---------------------------------------------------------
    def residuals(self, al, be, ga):
        """al,be,ga arrays (unwrapped, be_j<al_j<ga_j<be_{j+1}); returns
        relative residuals of the 3N conditions."""
        N = self.N
        res = []
        for j in range(N):
            tgt = np.log(self.R / self.t[j])
            res.append(self.int_side(al[j], be[j], al, be, ga) / tgt - 1)
            res.append(self.int_side(al[j], ga[j], al, be, ga) / tgt - 1)
        for j in range(N):
            jn = (j + 1) % N
            bnext = be[jn] + (TWO_PI if jn == 0 else 0.0)
            dth = (self.th[jn] - self.th[j]) % TWO_PI
            if N == 1:
                dth = TWO_PI
            res.append(self.int_arc(ga[j], bnext, al, be, ga) / dth - 1)
        return np.array(res)

    def logc(self, al, be, ga):
        s = 0.0
        for j in range(self.N):
            wf = lambda phi, a=al[j]: np.abs(phi - a)
            s += self.int_side(al[j], be[j], al, be, ga, wf)
            s += self.int_side(al[j], ga[j], al, be, ga, wf)
        return np.log(self.R) - s / TWO_PI


def solve_general(R, thetas, tips, n=64, npan=1, x0=None, verbose=0, steps=12):
    """Solve the parameter problem with no symmetry assumption.
    Unknowns: the 3N gaps between consecutive prevertices (positive, sum 2pi)
    -> parametrised by 3N logits with softmax; rotation fixed separately by an
    extra unknown offset.  Continuation in slit length from short slits."""
    thetas = np.asarray(thetas, float)
    tips = np.asarray(tips, float)
    N = len(thetas)
    Ltar = np.log(R / tips)

    def unpack(x):
        off = x[0]
        g = np.exp(x[1:] - np.max(x[1:]))
        g = g / g.sum() * TWO_PI
        pts = off + np.concatenate([[0.0], np.cumsum(g)[:-1]])
        be, al, ga = pts[0::3], pts[1::3], pts[2::3]
        return al, be, ga

    def make(lam):
        sd = SlitDisc(R, thetas, R * np.exp(-lam * Ltar), n=n, npan=npan)

        def fun(x):
            al, be, ga = unpack(x)
            r = sd.residuals(al, be, ga)
            # rotation fixing: prevertex of first tip relative to theta_0 is
            # not determined by lengths; pin sum of (alpha_j - theta_j) mean
            # via circular mean -> one extra equation
            rot = np.angle(np.mean(np.exp(1j * (al - thetas))))
            return np.concatenate([r, [rot]])
        return sd, fun

    if x0 is None:
        lam0 = 0.02
        # short slits: half-widths ~ slit length (half-plane slit of height l
        # has prevertices at +-l)
        hw = lam0 * Ltar
        pts = np.empty(3 * N)
        pts[0::3] = thetas - hw
        pts[1::3] = thetas
        pts[2::3] = thetas + hw
        gaps = np.diff(np.concatenate([pts, [pts[0] + TWO_PI]]))
        assert np.all(gaps > 0)
        x = np.concatenate([[pts[0]], np.log(gaps)])
        lams = np.concatenate([[lam0], np.linspace(lam0, 1.0, steps + 1)[1:]])
    else:
        x = x0
        lams = [1.0]
    for lam in lams:
        sd, fun = make(lam)
        sol = least_squares(fun, x, xtol=1e-15, ftol=1e-15, gtol=1e-15,
                            x_scale='jac', method='trf')
        x = sol.x
        if verbose:
            print("  lam=%.3f  max|res|=%.2e  nfev=%d" %
                  (lam, np.max(np.abs(sol.fun)), sol.nfev))
    al, be, ga = unpack(x)
    return sd, x, (al, be, ga), np.max(np.abs(sol.fun))


def conformal_radius(R, thetas, tips, n=64, npan=1, x0=None, verbose=0, steps=12):
    sd, x, (al, be, ga), rmax = solve_general(R, thetas, tips, n, npan, x0, verbose, steps)
    return np.exp(sd.logc(al, be, ga)), rmax, x, (al, be, ga)
