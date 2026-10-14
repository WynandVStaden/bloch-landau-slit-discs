"""Univalent Bloch-Landau constant: the classical explicit domains.

U is the infimum of inradius / conformal radius over simply connected domains,
so any domain with inradius <= 1 and conformal radius r at some point gives
U <= 1/r.

Domain here: the disc |z| < R with three radial slits [1, R] at angles 0, 120,
240 degrees and three radial slits [t, R] at 60, 180, 300 degrees. Its
conformal radius at 0 is explicit:

  * the unit disc minus the slit [s, 1] has conformal radius 4s/(1+s)^2 at 0,
    by composing the Koebe function k(z) = z/(1-z)^2 with a Moebius map;
  * z -> z^3 turns three symmetric slits into one;
  * the first map commutes with the reflections in the 60-degree lines, so the
    second set of slits is again radial after the first set is mapped away.
"""
import numpy as np
from scipy.optimize import minimize_scalar


def koebe(z):
    return z / (1 - z) ** 2


def slit_map_on_negative_axis(x, s):
    """Image of -x (0 < x < 1) under the conformal map of D minus [s,1] onto D fixing 0; returns beta with image -beta."""
    ks = koebe(s)
    w = koebe(-x)
    zeta = w / (1 - w / ks)
    kk = abs(zeta * (4 * ks + 1) / (4 * ks))  # k(-beta) = -kk, kk in (0, 1/4)
    return ((1 - 2 * kk) - np.sqrt(1 - 4 * kk)) / (2 * kk)


def conformal_radius(big_r, t):
    """Disc of radius big_r, slits [1, big_r] at 0/120/240 deg and [t, big_r] at 60/180/300 deg."""
    s3 = (1 / big_r) ** 3
    first = (4 * s3 / (1 + s3) ** 2) ** (1 / 3)
    beta = slit_map_on_negative_axis((t / big_r) ** 3, s3)
    second = (4 * beta / (1 + beta) ** 2) ** (1 / 3)
    return big_r * first * second


def inradius(big_r, t, n=600):
    """Numerical inradius (largest disc avoiding the slits and inside the big disc), by grid search plus refinement."""
    from scipy.optimize import minimize
    tips = []
    segs = []
    for k in range(3):
        a = 2 * np.pi * k / 3
        segs.append((np.exp(1j * a), big_r * np.exp(1j * a)))
        b = a + np.pi / 3
        segs.append((t * np.exp(1j * b), big_r * np.exp(1j * b)))

    def dist(p):
        z = complex(p[0], p[1])
        d = big_r - abs(z)
        for u, v in segs:
            w = v - u
            lam = np.clip(((z - u) * np.conj(w)).real / abs(w) ** 2, 0, 1)
            d = min(d, abs(z - (u + lam * w)))
        return d

    best = 0
    xs = np.linspace(0, big_r, n)
    # by symmetry it is enough to look in the sector 0..60 degrees
    for r in np.linspace(0, big_r, 160):
        for a in np.linspace(0, np.pi / 3, 60):
            z = r * np.exp(1j * a)
            d = dist((z.real, z.imag))
            if d > best - 0.05:
                res = minimize(lambda p: -dist(p), (z.real, z.imag), method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-12})
                best = max(best, -res.fun)
    return best


if __name__ == "__main__":
    for big_r, t in ((3.0, 2.0),):
        r = conformal_radius(big_r, t)
        print(f"R={big_r}, t={t}: inradius {inradius(big_r, t):.6f}, conformal radius {r:.7f}, bound U <= {1 / r:.7f}")
    print("known values: Robinson 0.658, Goodman 0.65647, Beller-Hummel 0.6564155, Carroll-Ortega-Cerda 0.6563937")
