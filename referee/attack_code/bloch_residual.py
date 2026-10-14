"""A posteriori check of the conformal-radius computation.

The solver represents the Green function as a sum of logarithmic potentials of
densities on the slits. Such a function is harmonic in the domain, has the
right singularity at 0 and vanishes on the circle automatically. So the only
thing that can be wrong is its value on the slits, which should be 0. This
script evaluates it at many points of every slit that were NOT used as
collocation points. By the maximum principle the error in log(conformal radius)
is at most the largest such value.

Usage: python bloch_residual.py R t2 angle3 t3 [angle4 t4 ...]
"""
import sys

import numpy as np

import bloch_design as d
from bloch_solver import phi


def solve(big_r, slits, n):
    big_l = np.log(big_r)
    k = np.arange(1, n + 1)
    x = np.cos((2 * k - 1) * np.pi / (2 * n))
    cheb = np.cos(np.outer(np.arccos(x), np.arange(n)))
    to_coef = (2 / n) * cheb.T
    to_coef[0] /= 2
    lam = np.pi * np.concatenate([[np.log(2)], 1 / np.arange(1, n)])
    pts = [s.w(x, big_l) for s in slits]

    def rows(i, xe):
        """Matrix rows giving the potential at parameters xe of slit i from the nodal values of all densities."""
        we = slits[i].w(xe, big_l)
        che = np.cos(np.outer(np.arccos(xe), np.arange(n)))
        blocks = []
        for j, sj in enumerate(slits):
            diff = we[:, None] - pts[j][None, :]
            if i == j:
                dx = xe[:, None] - x[None, :]
                same = np.abs(dx) < 1e-14
                with np.errstate(divide="ignore", invalid="ignore"):
                    smooth = phi(diff) + np.log(np.abs(dx))
                smooth[same] = np.broadcast_to(-np.log(np.abs(sj.dw(x)))[None, :], smooth.shape)[same]
                blocks.append((che * lam) @ to_coef + smooth * (np.pi / n))
            else:
                blocks.append(phi(diff) * (np.pi / n))
        return np.hstack(blocks)

    a_mat = np.vstack([rows(i, x) for i in range(len(slits))])
    rhs = np.concatenate([s.a * x for s in slits])
    psi = np.linalg.solve(a_mat, rhs)
    c1 = [(to_coef @ psi[j * n:(j + 1) * n])[1] for j in range(len(slits))]
    log_rho = big_l - (np.pi / 4) * sum(s.a * c for s, c in zip(slits, c1))
    return psi, rows, log_rho


def main():
    p = [float(v) for v in sys.argv[1:]]
    big_r = p[0]
    classes = [(0.0, 1.0), (np.pi / 3, p[1])] + [(p[2 + 2 * k], p[3 + 2 * k]) for k in range((len(p) - 2) // 2)]
    slits = d.slit_list(big_r, classes)
    xe = np.cos(np.linspace(0.003, np.pi - 0.003, 1501))  # dense test points, clustered at the tips
    for n in (40, 64, 96):
        psi, rows, log_rho = solve(big_r, slits, n)
        worst = 0.0
        for i, s in enumerate(slits):
            # u on the slit = (L - Re w) + potential = -a x + potential
            u = -s.a * xe + rows(i, xe) @ psi
            worst = max(worst, float(np.abs(u).max()))
        print(f"n={n:3d}: conformal radius {np.exp(log_rho):.13f}; largest |Green function| on the slits at 1501 test points each: {worst:.2e}")


if __name__ == "__main__":
    main()
