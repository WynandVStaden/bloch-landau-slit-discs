"""Floating-point preview of the conformal-radius certificate (not rigorous).

Freezes an odd, class-symmetric polynomial density from the solver, then evaluates the
resulting function U_h on each class representative with accurate quadrature
(many Gauss-Chebyshev nodes for the cross terms) at Chebyshev nodes, and reports
max |U_h| and the sum of absolute Chebyshev coefficients of the interpolant.

Usage: python prototype_residual.py [design] [n] [m] [N]
"""
import sys

import numpy as np

sys.path.insert(0, "..")
import bloch_curved2 as b  # noqa: E402
from bloch_residual import solve  # noqa: E402
from bloch_solver import phi  # noqa: E402
from designs import DESIGNS  # noqa: E402


def frozen_density(big_r, classes, n):
    """Odd Chebyshev coefficients per slit, averaged over each symmetry class. Returns slits, coef (len(slits) x n), class index."""
    slits = b.slit_objects(big_r, classes)
    psi, _, log_rho = solve(big_r, slits, n)
    k = np.arange(1, n + 1)
    x = np.cos((2 * k - 1) * np.pi / (2 * n))
    cheb = np.cos(np.outer(np.arccos(x), np.arange(n)))
    to_coef = (2 / n) * cheb.T
    to_coef[0] /= 2
    coef = np.array([to_coef @ psi[j * n:(j + 1) * n] for j in range(len(slits))])
    key = [(round(s.a, 12), round(abs(s.b1), 12), round(abs(s.b2), 12)) for s in slits]
    cls = [sorted(set(key)).index(v) for v in key]
    out = np.zeros_like(coef)
    for c in set(cls):
        members = [j for j in range(len(slits)) if cls[j] == c]
        mean = coef[members].mean(axis=0)
        mean[0::2] = 0.0  # odd part only
        out[members] = mean
    return slits, out, np.array(cls), log_rho


def evaluate(big_r, slits, coef, i, xe, n_quad):
    """U_h at parameters xe of slit i, from the frozen coefficients."""
    big_l = np.log(big_r)
    n = coef.shape[1]
    y = np.cos((2 * np.arange(1, n_quad + 1) - 1) * np.pi / (2 * n_quad))
    ty = np.cos(np.outer(np.arccos(y), np.arange(n)))
    lam = np.pi * np.concatenate([[np.log(2)], 1 / np.arange(1, n)])
    we = slits[i].w(xe, big_l)
    u = -slits[i].a * xe
    for j, sj in enumerate(slits):
        psi_y = ty @ coef[j]
        diff = we[:, None] - sj.w(y, big_l)[None, :]
        if i == j:
            dx = xe[:, None] - y[None, :]
            with np.errstate(divide="ignore", invalid="ignore"):
                smooth = phi(diff) + np.log(np.abs(dx))
            bad = np.abs(dx) < 1e-13
            smooth[bad] = np.broadcast_to(-np.log(np.abs(sj.dw(y)))[None, :], smooth.shape)[bad]
            u += np.cos(np.outer(np.arccos(xe), np.arange(n))) @ (lam * coef[j]) + smooth @ psi_y * (np.pi / n_quad)
        else:
            u += phi(diff) @ psi_y * (np.pi / n_quad)
    return u


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "G"
    design = next(v for k, v in DESIGNS.items() if k.startswith(name))
    m = int(sys.argv[3]) if len(sys.argv) > 3 else 1500
    n_quad = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
    for n in ([int(sys.argv[2])] if len(sys.argv) > 2 and sys.argv[2] != "-" else [32, 48, 64, 96, 128]):
        slits, coef, cls, log_rho = frozen_density(design["R"], design["classes"], n)
        a = np.array([s.a for s in slits])
        log_rho_frozen = np.log(design["R"]) - (np.pi / 4) * float(a @ coef[:, 1])
        xe = np.cos((2 * np.arange(1, m + 2) - 1) * np.pi / (2 * (m + 1)))
        line = []
        for c in sorted(set(cls)):
            i = int(np.argmax(cls == c))
            u = evaluate(design["R"], slits, coef, i, xe, n_quad)
            cc = (2 / (m + 1)) * np.cos(np.outer(np.arange(m + 1), np.arccos(xe))) @ u
            cc[0] /= 2
            line.append(f"[a={slits[i].a:.3f}: max {np.abs(u).max():.1e}, sum|c| {np.abs(cc).sum():.1e}, |psi| {np.abs(coef[i]).sum():.2f}]")
        print(f"n={n:3d}: log rho solver {log_rho:.13f}, frozen {log_rho_frozen:.13f}, rho {np.exp(log_rho_frozen):.12f}\n   " + "\n   ".join(line), flush=True)
