import sys, time
import numpy as np
from sc_slitmap import conformal_radius
from geom import R, DEG, claimed_slits

np.set_printoptions(precision=12, linewidth=150)


def calib1():
    print("== Calibration 1: unit disc minus [s,1]; exact 4s/(1+s)^2")
    for s in (0.5, 0.2, 0.05):
        for n, npan in ((32, 1), (64, 1), (64, 4)):
            c, rmax, _, _ = conformal_radius(1.0, [0.0], [s], n=n, npan=npan)
            ex = 4 * s / (1 + s) ** 2
            print("  s=%.2f n=%d npan=%d  c=%.15f exact=%.15f  err=%.2e  res=%.1e" %
                  (s, n, npan, c, ex, c - ex, rmax))


def calib2():
    print("== Calibration 2: R=3, [1,3] at 0,120,240 and [2,3] at 60,180,300; expected 1.5207430363")
    th = np.array([0, 60, 120, 180, 240, 300]) * DEG
    tp = [1, 2, 1, 2, 1, 2]
    L = 27 + 1 / 27 + 27 / 8 + 8 / 27
    ex = (108 / L) ** (1 / 3)
    for n, npan in ((32, 1), (64, 1), (64, 4), (100, 6)):
        c, rmax, _, _ = conformal_radius(3.0, th, tp, n=n, npan=npan)
        print("  n=%d npan=%d  c=%.15f  exact(Joukowski)=%.15f err=%.2e res=%.1e" %
              (n, npan, c, ex, c - ex, rmax))


def claimed():
    print("== Claimed 18-slit domain, R=%.5f (no symmetry imposed: 54 free prevertices)" % R)
    sl = sorted(claimed_slits())
    th = np.array([a for a, t in sl])
    tp = np.array([t for a, t in sl])
    out = []
    x0 = None
    for n, npan in ((40, 1), (64, 1), (64, 3), (96, 5), (128, 8)):
        t0 = time.time()
        c, rmax, x0, (al, be, ga) = conformal_radius(R, th, tp, n=n, npan=npan, x0=x0,
                                                     verbose=(x0 is None))
        out.append(c)
        print("  n=%d npan=%d  crad=%.15f  max rel residual=%.1e  1/crad=%.12f  (%.1fs)" %
              (n, npan, c, rmax, 1 / c, time.time() - t0))
    print("  claimed crad 1.523513988544 ; diff of finest = %.3e" % (out[-1] - 1.523513988544))
    print("  spread among last three resolutions: %.2e" % (max(out[-3:]) - min(out[-3:])))
    print("  ratio 1/crad = %.12f   (C-OC bound 0.6563937)" % (1 / out[-1]))
    print("  harmonic measure at 0 of outer-circle arcs: %.6f" %
          (np.sum(np.concatenate([be[1:], [be[0] + 2 * np.pi]]) - ga) / (2 * np.pi)))
    print("  prevertices (deg) be, al, ga for first 6 slits:")
    for j in range(6):
        print("   slit %.5f tip %.5f: %.8f %.8f %.8f" %
              (th[j] / DEG, tp[j], be[j] / DEG, al[j] / DEG, ga[j] / DEG))


if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3"]
    if "1" in which: calib1()
    if "2" in which: calib2()
    if "3" in which: claimed()
