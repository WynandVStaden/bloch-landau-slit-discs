"""Rigorous lower bound for the conformal radius of a slit disc.

THE ARGUMENT
------------
Work in w = log z, L = log R. Slit j together with its reflection in Re w = L is the arc
    w_j(x) = L + a_j x + i (theta_j + b1_j (a_j x)^2 + b2_j (a_j x)^4),   -1 <= x <= 1,
and the slit itself is -1 <= x <= 0. Let psi_j be ODD polynomials, the same for
slits that are images of one another under the symmetry group, with rational
Chebyshev coefficients c_{j,n} (they come from the floating-point solver, but
nothing below assumes they solve anything). Put
    U(w) = (L - Re w) + sum_j  int_{-1}^{1} Phi(w - w_j(y)) psi_j(y) dy / sqrt(1 - y^2),
    Phi(w) = -log|2 sinh(w/2)|.
Then, exactly:
  * U, as a function of z = e^w, is harmonic in D minus {0} and continuous up to the boundary;
  * U = 0 on |z| = R (oddness of psi_j);
  * U = -log|z| + c + o(1) at z = 0 with c = L - (pi/4) sum_j a_j c_{j,1}.
If U <= delta on every slit, the maximum principle applied to g - (U - delta),
g the Green function of D, gives
    log(conformal radius) >= c - delta.

BOUNDING U ON A SLIT
--------------------
f_i(x) = U(w_i(x)) is odd in x. We bound sup |f_i| over [-1, 1] by
    sum of |Chebyshev coefficients of the interpolant p of f_i at M Chebyshev nodes|
    + 2 * sum_{k >= M} |a_k|,   (aliasing bound for interpolation at the zeros of T_M)
where a_k are the Chebyshev coefficients of f_i. f_i is the real part, on
[-1, 1], of a function analytic in a neighbourhood of a Bernstein ellipse
E_rho, so |a_k| <= 4 B / (rho^k + rho^-k) with B >= sup over E_rho of the real
part. That real part is made of explicit polynomials and of integrals of
log|kernel|, and |kernel| is bounded above and below, away from 0, on a
rectangle containing E_rho by a branch-and-bound in ball arithmetic.

The values f_i(x_k) are computed in ball arithmetic. The logarithmic part of
the self-interaction is integrated exactly; every other integral is a
Gauss-Chebyshev sum plus a remainder bounded in the same way (analyticity of
the kernel in the integration variable on a rectangle, verified by branch-and-bound).

All trusted computations use Arb (python-flint). Floating point is used only to
choose parameters (ellipse sizes, numbers of nodes), to guide the branch and bound, and to
produce the c_{j,n}. Enclosures are passed between processes as a double-precision midpoint
and an upward-rounded radius, which is exact.

Usage: python conformal_certificate.py [design, default G] [n, default 64] [processes, default 8]
"""
import json
import math
import os
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

import numpy as np
from flint import acb, arb, ctx, fmpq

from inradius_certificate import DESIGNS, PI3, arb_angle, arb_frac, float_angle, slit_table

DESIGNS = dict(DESIGNS)
DESIGNS["R6"] = ("3", [("0", "1", "0", "0"), (PI3, "2", "0", "0")])  # exactly solvable test case
DESIGNS["ONE"] = ("1", [("0", "0.5", "0", "0")])  # placeholder, replaced below: one slit only

PREC = 100
TOL_INTERP = 1e-9
TOL_QUAD = 1e-13
SAFETY = 0.5  # fraction of the estimated distance to the nearest singularity used for the rectangles
I = acb(0, 1)


def up(x):
    """A float that is certainly >= the Arb number x (x must be non-negative)."""
    v = float(x.upper())
    assert v >= 0
    return v * (1 + 2.0**-40) + 1e-300


class Geometry:
    """The arcs in ball arithmetic."""

    def __init__(self, name):
        ctx.prec = PREC
        r_str, classes = DESIGNS[name]
        self.big_r = Fraction(r_str)
        table = slit_table(classes)
        if name == "ONE":
            table = table[:1]
        self.big_l = arb_frac(self.big_r).log()
        self.cls, self.a, self.th, self.b1, self.b2 = [], [], [], [], []
        tips = [Fraction(c[1]) for c in classes]
        for r, sgn, alpha, tip, b1, b2 in table:
            self.cls.append(tips.index(tip))
            self.a.append((arb_frac(self.big_r) / arb_frac(tip)).log())
            self.th.append(arb_angle(r, sgn, alpha))
            self.b1.append(arb_frac(b1))
            self.b2.append(arb_frac(b2))
        self.n_slits = len(table)
        self.table = table
        # floating-point copies, for planning only
        self.fa = [float(v.mid()) for v in self.a]
        self.fth = [float_angle(r, sgn, alpha) for r, sgn, alpha, *_ in table]
        self.fb1 = [float(t[4]) for t in table]
        self.fb2 = [float(t[5]) for t in table]
        # sup |w_j'| on a rectangle |Re x| <= 1.2, |Im x| <= 0.6 (crude, for scaling the branch-and-bound only)
        self.speed = [a * (1 + 2 * abs(b1) * a * 1.4 + 4 * abs(b2) * a**3 * 1.4**3)
                      for a, b1, b2 in zip(self.fa, self.fb1, self.fb2)]
        self.reps = [self.cls.index(c) for c in sorted(set(self.cls))]

    def w(self, j, x):
        s = self.a[j] * x
        s2 = s * s
        return self.big_l + s + I * (self.th[j] + self.b1[j] * s2 + self.b2[j] * s2 * s2)

    def cross(self, i, j, x, y):
        """2 sinh((w_i(x) - w_j(y)) / 2)."""
        return 2 * ((self.w(i, x) - self.w(j, y)) / 2).sinh()

    def self_kernel(self, i, x, y):
        """2 sinh((w_i(x) - w_i(y)) / 2) / (x - y), written without the removable singularity."""
        a = self.a[i]
        dd = a + I * (self.b1[i] * a**2 * (x + y) + self.b2[i] * a**4 * (x * x * x + x * x * y + x * y * y + y * y * y))
        u = (x - y) * dd
        return dd * (I * u / 2).sinc()  # sinh(t)/t = sinc(i t)

    def fw(self, j, x):
        s = self.fa[j] * x
        return math.log(float(self.big_r)) + s + 1j * (self.fth[j] + self.fb1[j] * s**2 + self.fb2[j] * s**4)


def bound_abs(func, box, scale, ratio=1e3, min_width=1e-7, max_boxes=400000):
    """Rigorous (lower, upper) bounds, as Arb numbers, of |func(x, y)| for x = xr + i xi, y = yr + i yi in the box.

    box: four (lo, hi) pairs of floats. Returns None if 0 cannot be excluded.
    """
    lo_all, hi_all = None, None
    stack, count = [box], 0
    while stack:
        bx = stack.pop()
        count += 1
        if count > max_boxes:
            return None
        balls = [arb(lo).union(arb(hi)) for lo, hi in bx]  # exact endpoints, so the ball contains the interval
        v = func(acb(balls[0], balls[1]), acb(balls[2], balls[3]))
        lo, hi = v.abs_lower(), v.abs_upper()
        if lo > 0 and hi < ratio * lo:
            lo_all = lo if lo_all is None or lo < lo_all else lo_all
            hi_all = hi if hi_all is None or hi > hi_all else hi_all
            continue
        widths = [(h - l) * s for (l, h), s in zip(bx, scale)]
        d = int(np.argmax(widths))
        if widths[d] < min_width:
            return None
        mid = (bx[d][0] + bx[d][1]) / 2
        for part in ((bx[d][0], mid), (mid, bx[d][1])):
            child = list(bx)
            child[d] = part
            stack.append(child)
    return lo_all, hi_all


def rho_of(b):
    return b + math.sqrt(1 + b * b)


def make_density(geo, n):
    """Untrusted floating-point solve; returns odd Chebyshev coefficients per class as exact doubles."""
    sys.path.insert(0, "..")
    import bloch_curved2 as bc
    from bloch_residual import solve

    slits = [bc.Slit2(geo.fth[j], geo.fa[j], geo.fb1[j], geo.fb2[j]) for j in range(geo.n_slits)]
    psi, _, log_rho = solve(float(geo.big_r), slits, n)
    k = np.arange(1, n + 1)
    x = np.cos((2 * k - 1) * np.pi / (2 * n))
    to_coef = (2 / n) * np.cos(np.outer(np.arccos(x), np.arange(n))).T
    to_coef[0] /= 2
    coef = np.array([to_coef @ psi[j * n:(j + 1) * n] for j in range(geo.n_slits)])
    out = {}
    for c in sorted(set(geo.cls)):
        mean = coef[[j for j in range(geo.n_slits) if geo.cls[j] == c]].mean(axis=0)
        mean[0::2] = 0.0
        out[c] = [float(v) for v in mean]
    return out, float(np.exp(log_rho))


def float_distance(geo, i, j, x_range):
    """Estimated distance on the cylinder between arc i (parameters in x_range) and the whole arc j."""
    xs = np.linspace(x_range[0], x_range[1], 600)
    ys = np.linspace(-1, 1, 600)
    wi = np.array([geo.fw(i, v) for v in xs])
    wj = np.array([geo.fw(j, v) for v in ys])
    d = wi[:, None] - wj[None, :]
    im = (d.imag + np.pi) % (2 * np.pi) - np.pi
    return float(np.hypot(d.real, im).min())


def ellipse(b):
    """For a rectangle half-height b (a double): a half-width a_half and an ellipse parameter r (doubles) with
    E_r contained in [-a_half, a_half] x [-b, b]. The inclusion is verified in ball arithmetic."""
    a_half = math.sqrt(1 + b * b) * (1 + 1e-12)
    r = rho_of(b) * (1 - 1e-12)
    ra = arb(r)
    assert ra > 1 and (ra + 1 / ra) / 2 < arb(a_half) and (ra - 1 / ra) / 2 < arb(b), "ellipse not inside the rectangle"
    return a_half, r


def log_sup(bounds):
    """An Arb upper bound for sup |log |F|| from Arb (lower, upper) bounds of |F|."""
    lo, hi = bounds
    return abs(lo.log()).max(abs(hi.log()))


def plan(geo, coef, log=print):
    """Choose rectangles and node counts, and prove the kernel bounds they need.

    Floating point is used only to CHOOSE b, the node counts M and N and the boxes of the branch and bound.
    Every bound (psi_sup, A, B, the interpolation remainder, the quadrature remainders) is computed in ball
    arithmetic; it is converted to a double, rounded upwards, only when it is handed to the worker processes.
    """
    for vals in coef.values():
        assert all(v == 0.0 for v in vals[0::2]), "the density polynomials must be odd"
    deg = max(n for c in coef.values() for n, v in enumerate(c) if v != 0.0)
    psi_sup = {c: sum((abs(arb(v)) for v in vals), arb(0)) for c, vals in coef.items()}  # >= sup |psi| on [-1, 1]
    pi = arb.pi()
    out = {}
    for i in geo.reps:
        t0 = time.time()
        others = [j for j in range(geo.n_slits) if j != i]
        # ---- interpolation in x: rectangle around [-1, 1] for arc i
        dist = min([float_distance(geo, i, j, (-1, 1)) for j in others] + [9.0])
        b = min(0.3, SAFETY * dist / geo.speed[i])
        while True:
            a_half, r = ellipse(b)
            box = [(-a_half, a_half), (-b, b), (-1.0, 1.0), (0.0, 0.0)]
            scale = [geo.speed[i], geo.speed[i], 1.0, 0.0]
            big_b, ok = arb(0), True
            for j in [i] + others:
                func = (lambda x, y, j=j: geo.self_kernel(i, x, y)) if j == i else (lambda x, y, j=j: geo.cross(i, j, x, y))
                scale[2] = geo.speed[j]
                bd = bound_abs(func, box, scale)
                if bd is None:
                    ok = False
                    break
                big_b += pi * psi_sup[geo.cls[j]] * log_sup(bd)
            if ok:
                break
            b /= 2
        ra = arb(r)
        big_b += geo.a[i] * (ra + 1 / ra) / 2
        for n, v in enumerate(coef[geo.cls[i]]):
            if v != 0.0:
                big_b += abs(arb(v)) * (pi / n) * (ra**n + ra ** (-n)) / 2

        def interp_remainder(m):
            return 8 * big_b * ra ** (1 - m) / (ra - 1)

        m_nodes = max(2, 2 * int(math.log(8 * up(big_b) / ((r - 1) * TOL_INTERP)) / math.log(r) / 2))  # floating-point guess
        while not interp_remainder(m_nodes) < TOL_INTERP:
            m_nodes += 2
        interp_err = up(interp_remainder(m_nodes))
        # ---- quadrature in y: rectangle around [-1, 1] for each source arc j, target points on the slit x in [-1, 0]
        quad = {}
        for j in [i] + others:
            if j == i:
                bq = 0.3
            else:
                bq = min(0.3, SAFETY * float_distance(geo, i, j, (-1, 0)) / geo.speed[j])
            while True:
                a_half_q, rq = ellipse(bq)
                box = [(-1.0, 0.0), (0.0, 0.0), (-a_half_q, a_half_q), (-bq, bq)]
                func = (lambda x, y: geo.self_kernel(i, x, y)) if j == i else (lambda x, y, j=j: geo.cross(i, j, x, y))
                bd = bound_abs(func, box, [geo.speed[i], 0.0, geo.speed[j], geo.speed[j]])
                if bd is not None:
                    break
                bq /= 2
            rqa = arb(rq)
            const = 8 * pi * psi_sup[geo.cls[j]] * log_sup(bd) / (rqa - 1)

            def quad_remainder(n_q):
                return const * rqa ** (-(2 * n_q - 1 - deg))

            n_q = deg // 2 + 2
            n_q = max(n_q, int((deg + 1 + math.log(max(up(const), 1e-300) / TOL_QUAD) / math.log(rq)) / 2))  # guess
            while not quad_remainder(n_q) < TOL_QUAD:
                n_q += 1
            quad[j] = (n_q, up(quad_remainder(n_q)))
        out[i] = {"b": b, "rho": r, "B": up(big_b), "M": m_nodes, "interp_err": interp_err, "quad": quad}
        log(f"  slit {i} (class {geo.cls[i]}, a = {geo.fa[i]:.4f}): r = {r:.5f}, B = {up(big_b):.1f}, {m_nodes} interpolation nodes, "
            f"quadrature nodes {min(q[0] for q in quad.values())}..{max(q[0] for q in quad.values())}, "
            f"{m_nodes // 2 * sum(q[0] for q in quad.values()):.2e} kernel evaluations; planned in {time.time() - t0:.0f} s")
    return out


# ------------------------------------------------------------------ evaluation of f_i at the nodes (worker processes)

_STATE = {}


def _init(name, coef):
    ctx.prec = PREC
    _STATE["geo"] = Geometry(name)
    _STATE["coef"] = {int(c): v for c, v in coef.items()}
    _STATE["nodes"] = {}


def _source_nodes(j, n_q):
    """w_j(y_l) and psi_j(y_l) * pi / n_q at the n_q Gauss-Chebyshev nodes, cached."""
    geo, coef = _STATE["geo"], _STATE["coef"][_STATE["geo"].cls[j]]
    key = (j, n_q)
    if key not in _STATE["nodes"]:
        cos_tab = [arb.cos_pi_fmpq(fmpq(t, 2 * n_q)) for t in range(4 * n_q)]
        ws, ps, ys = [], [], []
        wt = arb.pi() / n_q
        for l in range(1, n_q + 1):
            y = cos_tab[2 * l - 1]
            psi = arb(0)
            for n, c in enumerate(coef):
                if c != 0.0:
                    psi += arb(c) * cos_tab[(n * (2 * l - 1)) % (4 * n_q)]
            ws.append(geo.w(j, y)), ps.append(psi * wt), ys.append(y)
        _STATE["nodes"][key] = (ws, ps, ys)
    return _STATE["nodes"][key]


def _eval_nodes(task):
    i, m_nodes, ks, quad = task
    geo, coef = _STATE["geo"], _STATE["coef"][_STATE["geo"].cls[i]]
    out = []
    for k in ks:
        x = arb.cos_pi_fmpq(fmpq(2 * k + 1, 2 * m_nodes))
        wx = geo.w(i, x)
        f = -geo.a[i] * x
        # exact logarithmic part of the self-interaction: sum_n c_n (pi / n) T_n(x)
        for n, c in enumerate(coef):
            if c != 0.0:
                f += arb(c) * (arb.pi() / n) * arb.cos_pi_fmpq(fmpq(n * (2 * k + 1), 2 * m_nodes))
        for j, (n_q, err) in quad.items():
            ws, ps, ys = _source_nodes(j, n_q)
            acc = arb(0)
            if j == i:
                for y, p in zip(ys, ps):
                    acc -= abs(geo.self_kernel(i, x, y)).log() * p
            else:
                for wy, p in zip(ws, ps):
                    acc -= abs(((wx - wy) / 2).sinh() * 2).log() * p
            f += acc + arb(0, err)
        mid = float(f.mid())
        rad = up(f.rad()) + abs(mid) * 2.0**-50
        out.append((k, mid, rad))
    return i, out


def certify(name, n, processes, log=print):
    ctx.prec = PREC
    t_start = time.time()
    geo = Geometry(name)
    # The density is frozen in a data file so that the proof does not depend on the floating-point solver:
    # the argument is valid for ANY odd coefficients; they only influence how small delta comes out.
    density_file = f"density_{name}_n{n}.json"
    if os.path.exists(density_file):
        with open(density_file) as fh:
            stored = json.load(fh)
        coef, rho_float = {int(c): v for c, v in stored["coefficients"].items()}, stored["floating_point_conformal_radius"]
        log(f"density coefficients read from {density_file}")
    else:
        coef, rho_float = make_density(geo, n)
        with open(density_file, "w") as fh:
            json.dump({"design": name, "n": n, "coefficients": coef, "floating_point_conformal_radius": rho_float}, fh, indent=1)
        log(f"density coefficients computed by the floating-point solver and written to {density_file}")
    log(f"design {name}: {geo.n_slits} slits, solver size n = {n}; floating-point conformal radius {rho_float:.12f}")
    log("planning (branch-and-bound proofs of the kernel bounds):")
    pl = plan(geo, coef, log)
    delta = arb(0)
    with Pool(processes, initializer=_init, initargs=(name, coef)) as pool:
        for i in geo.reps:
            t0 = time.time()
            m_nodes = pl[i]["M"]
            ks = list(range(m_nodes // 2, m_nodes))  # nodes with x_k < 0: the slit itself
            chunk = max(1, len(ks) // (processes * 6))
            tasks = [(i, m_nodes, ks[s:s + chunk], pl[i]["quad"]) for s in range(0, len(ks), chunk)]
            vals = {}
            for _, res in pool.imap_unordered(_eval_nodes, tasks):
                for k, mid, rad in res:
                    vals[k] = arb(mid, rad)
            # interpolant of the odd function f_i at the M nodes: coefficients for odd q only
            cos_tab = [arb.cos_pi_fmpq(fmpq(t, 2 * m_nodes)) for t in range(4 * m_nodes)]
            total, worst = arb(0), arb(0)
            for q in range(1, m_nodes, 2):
                s = arb(0)
                for k in ks:
                    s += vals[k] * cos_tab[(q * (2 * k + 1)) % (4 * m_nodes)]
                total += abs(s) * 4 / m_nodes
            for k in ks:
                worst = worst.max(abs(vals[k]))
            d_i = total + arb(pl[i]["interp_err"])
            delta = delta.max(d_i)
            log(f"  slit {i} (class {geo.cls[i]}): max |f| at nodes <= {up(worst):.3e}, sum of |coefficients| <= {up(total):.3e}, "
                f"interpolation remainder <= {pl[i]['interp_err']:.1e}; sup |U| <= {up(d_i):.3e}   [{time.time() - t0:.0f} s]")
    c1 = arb(0)
    for j in range(geo.n_slits):
        c1 += geo.a[j] * arb(coef[geo.cls[j]][1])
    assert delta >= 0
    log_rho_lower = geo.big_l - arb.pi() / 4 * c1 - delta
    rho_lower = log_rho_lower.exp()
    lower = math.floor(float(rho_lower.lower()) * (1 - 2.0**-40) * 1e12) / 1e12  # rounded down to 12 decimals
    assert arb(lower) < rho_lower
    log('RHO_LOWER_BALL ' + rho_lower.lower().str(25) + ' DELTA ' + delta.upper().str(12))
    log(f"delta = {up(delta):.3e}")
    log(f"c = L - (pi/4) sum a_j c_j1 = {float((geo.big_l - arb.pi() / 4 * c1).mid()):.13f}")
    log(f"CERTIFIED: conformal radius of design {name} >= {lower:.12f}   (floating-point value {rho_float:.12f}); "
        f"total time {time.time() - t_start:.0f} s")
    with open(f"conformal_certificate_{name}.json", "w") as fh:
        json.dump({"design": name, "n": n, "coefficients": coef, "delta_upper": up(delta), "conformal_radius_lower": lower}, fh, indent=1)
    return lower


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "G"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    procs = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    lower = certify(name, n, procs)
    if name == "R6":
        sys.path.insert(0, "..")
        import bloch_explicit as be

        exact = be.conformal_radius(3.0, 2.0)
        print(f"exact value (closed form) {exact:.12f}; certified lower bound is below it by {exact - lower:.2e}")
    if name == "ONE":
        print(f"exact value 4s/(1+s)^2 = {4 * 0.5 / 1.5**2:.12f}; certified lower bound is below it by {4 * 0.5 / 1.5**2 - lower:.2e}")
