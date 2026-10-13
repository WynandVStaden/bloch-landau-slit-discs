"""Walk-on-spheres Monte Carlo sanity check of the conformal radius.
log crad = E_0[ log |B_tau| ]  (B = Brownian motion from 0, tau = exit time of D),
because g(z,0) + log|z| is harmonic with boundary values log|z|.
Uses the actual 2-D geometry (exact point-segment distances), no symmetry,
no conformal-map formula -> independent of sc_slitmap.py.  Low accuracy."""
import sys
import numpy as np
from geom import R, claimed_slits, dist_to_boundary


def wos(slits, R, nwalk, eps, rng):
    z = np.zeros(nwalk, complex)
    out = np.empty(nwalk)
    alive = np.arange(nwalk)
    it = 0
    while alive.size:
        d, _ = dist_to_boundary(z[alive], slits, R)
        done = d < eps
        out[alive[done]] = np.log(np.abs(z[alive[done]]))
        alive = alive[~done]
        d = d[~done]
        z[alive] += d * np.exp(2j * np.pi * rng.random(alive.size))
        it += 1
    return out


if __name__ == "__main__":
    rng = np.random.default_rng(20261013)
    nw = int(float(sys.argv[1])) if len(sys.argv) > 1 else 2_000_000
    # calibration 2
    cal = [(k * np.pi / 3, 1.0 if k % 2 == 0 else 2.0) for k in range(6)]
    for name, sl, RR, ref in (("calibration R=3 six slits", cal, 3.0, 1.5207430363),
                              ("claimed 18-slit domain", claimed_slits(), R, 1.523513988544)):
        for eps in (1e-4, 1e-6):
            tot = []
            for chunk in range(nw // 500_000):
                tot.append(wos(sl, RR, 500_000, eps, rng))
            v = np.concatenate(tot)
            m, se = v.mean(), v.std() / np.sqrt(v.size)
            print("%s eps=%.0e N=%d: crad = %.5f +- %.5f (1 sigma)   reference %.6f  z-score %.2f"
                  % (name, eps, v.size, np.exp(m), np.exp(m) * se, ref,
                     (np.exp(m) - ref) / (np.exp(m) * se)), flush=True)
