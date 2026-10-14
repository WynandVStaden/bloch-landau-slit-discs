"""Conformal radius of a disc with slits, to high accuracy.

Domain: the disc |z| < R minus slits that start on the circle |z| = R, leave it
at a right angle, and end at a free tip inside. Work in w = log z, where the
disc becomes the half-plane Re w < L = log R (2*pi*i-periodic) and a radial
slit becomes a horizontal segment.

The Green function g(z, 0), written as u(w), vanishes on Re w = L and on the
slits, so it extends to an odd function across Re w = L. Each slit together
with its mirror image is then one analytic arc with two free tips,

    w_j(x) = L + a_j x + i (alpha_j + c_j (a_j x)^2),   -1 <= x <= 1,

and u = (L - Re w) + sum_j  int  Phi(w - w_j(y)) psi_j(y) / sqrt(1 - y^2) dy,
Phi(w) = -log|2 sinh(w/2)|, with psi_j smooth. The condition u = 0 on the arcs
is a first-kind integral equation with a logarithmic kernel, solved by the
standard Chebyshev method (the log part is diagonal in the Chebyshev basis,
the remainder is smooth and integrated by Gauss-Chebyshev quadrature), which
converges exponentially. Finally

    log(conformal radius at 0) = L - (pi/4) * sum_j a_j * (first Chebyshev coefficient of psi_j).
"""
import numpy as np


class Slit:
    def __init__(self, alpha, depth, curve=0.0):
        """alpha: polar angle where the slit meets the circle; depth = log(R / tip radius); curve: bending coefficient."""
        self.alpha, self.a, self.c = alpha, depth, curve

    def w(self, x, big_l):
        s = self.a * x
        return big_l + s + 1j * (self.alpha + self.c * s * s)

    def dw(self, x):
        s = self.a * x
        return self.a * (1 + 2j * self.c * s)

    def z_points(self, big_r, n=400):
        """Points of the slit in the z-plane, from the circle to the tip."""
        x = -np.linspace(0, 1, n)
        return np.exp(self.w(x, np.log(big_r)))


def phi(w):
    """-log|2 sinh(w/2)|, the 2*pi*i-periodic fundamental solution."""
    return -np.log(np.abs(2 * np.sinh(w / 2)))


def conformal_radius(big_r, slits, n=48):
    big_l = np.log(big_r)
    k = np.arange(1, n + 1)
    x = np.cos((2 * k - 1) * np.pi / (2 * n))  # Chebyshev nodes, used for collocation and quadrature
    cheb = np.cos(np.outer(np.arccos(x), np.arange(n)))  # T_m(x_i)
    # psi values at nodes -> Chebyshev coefficients
    to_coef = (2 / n) * cheb.T
    to_coef[0] /= 2
    # log-kernel operator acting on coefficients: int -log|x-y| T_m(y)/sqrt(1-y^2) dy
    log_op = cheb * (np.pi * np.concatenate([[np.log(2)], 1 / np.arange(1, n)]))
    j_count = len(slits)
    a_mat = np.zeros((j_count * n, j_count * n))
    rhs = np.zeros(j_count * n)
    pts = [s.w(x, big_l) for s in slits]
    for i, si in enumerate(slits):
        rhs[i * n:(i + 1) * n] = si.a * x  # -(L - Re w) on the arc
        for j, sj in enumerate(slits):
            diff = pts[i][:, None] - pts[j][None, :]
            if i == j:
                with np.errstate(divide="ignore", invalid="ignore"):
                    smooth = phi(diff) + np.log(np.abs(x[:, None] - x[None, :]))
                np.fill_diagonal(smooth, -np.log(np.abs(sj.dw(x))))
                block = log_op @ to_coef + smooth * (np.pi / n)
            else:
                block = phi(diff) * (np.pi / n)
            a_mat[i * n:(i + 1) * n, j * n:(j + 1) * n] = block
    psi = np.linalg.solve(a_mat, rhs)
    c1 = np.array([(to_coef @ psi[j * n:(j + 1) * n])[1] for j in range(j_count)])
    return float(np.exp(big_l - (np.pi / 4) * sum(s.a * c for s, c in zip(slits, c1))))


def symmetric_slits(order, alpha, depth, curve=0.0, mirror=True):
    """Copies of a slit under rotation by 2*pi/order and, optionally, reflection in the real axis."""
    out = []
    for r in range(order):
        rot = 2 * np.pi * r / order
        out.append(Slit(rot + alpha, depth, curve))
        if mirror and abs(np.sin(order * alpha / 1.0)) > 1e-12 and abs(np.sin(order * alpha / 2.0)) > 1e-12:
            out.append(Slit(rot - alpha, depth, -curve))
    return out


if __name__ == "__main__":
    import bloch_explicit as be

    print("check 1: unit disc with one radial slit [s, 1]; exact conformal radius 4s/(1+s)^2")
    for s in (0.3, 0.6, 0.9):
        for n in (16, 32, 64):
            got = conformal_radius(1.0, [Slit(0.0, -np.log(s))], n)
            print(f"   s={s}  n={n:3d}  computed {got:.14f}  exact {4 * s / (1 + s) ** 2:.14f}  error {abs(got - 4 * s / (1 + s) ** 2):.1e}")
    print("check 2: disc of radius R with slits [1,R] at 0/120/240 and [t,R] at 60/180/300 degrees")
    for big_r, t in ((3.0, 2.0), (2.9562, 2.0)):
        exact = be.conformal_radius(big_r, t)
        slits = symmetric_slits(3, 0.0, np.log(big_r)) + symmetric_slits(3, np.pi / 3, np.log(big_r / t))
        for n in (24, 48, 96):
            got = conformal_radius(big_r, slits, n)
            print(f"   R={big_r} t={t}  n={n:3d}  computed {got:.13f}  exact {exact:.13f}  error {abs(got - exact):.1e}")
