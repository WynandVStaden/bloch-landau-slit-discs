# Referee: independent Lipschitz quadtree over the WHOLE disc (no symmetry used).
# For a square cell with centre c and half-diagonal rho:  sup_cell dist(.,bdry) <= |c-q| + rho for any boundary point q.
# Cells with |c-q| + rho <= THR are discarded. Double precision (errors ~1e-15 << margins); not a formal proof, an independent check.
import sys, numpy as np, time
from geom import *
from scipy.spatial import cKDTree
THR = 1 + float(sys.argv[1]) if len(sys.argv) > 1 else 1 + 1e-9
pts, sid, spar = [], [], []
for k, sl in enumerate(SL_F):
    s = np.linspace(-sl['a'], 0, 20001); pts.append(gam_f(sl, s)); sid.append(np.full(len(s), k)); spar.append(s)
P = np.concatenate(pts); SID = np.concatenate(sid); SPAR = np.concatenate(spar)
tree = cKDTree(np.c_[P.real, P.imag])
A = np.array([s['a'] for s in SL_F]); TH = np.array([s['th'] for s in SL_F]); B1 = np.array([s['b1'] for s in SL_F]); B2 = np.array([s['b2'] for s in SL_F])
def dist(c):
    d0, i = tree.query(np.c_[c.real, c.imag]); k = SID[i]; s = SPAR[i].copy()
    a, th, b1, b2 = A[k], TH[k], B1[k], B2[k]
    for _ in range(8):
        g = R_F*np.exp(s + 1j*(th + b1*s**2 + b2*s**4)); p1 = 1 + 1j*(2*b1*s + 4*b2*s**3)
        g1 = g*p1; g2 = g*(p1*p1 + 1j*(2*b1 + 12*b2*s*s))
        f = ((g - c)*np.conj(g1)).real; fp = np.abs(g1)**2 + ((g - c)*np.conj(g2)).real
        s = np.clip(np.where(fp > 0, s - f/np.where(fp > 0, fp, 1), s), -a, 0)
    g = R_F*np.exp(s + 1j*(th + b1*s**2 + b2*s**4))
    return np.minimum(np.minimum(d0, np.abs(g - c)), R_F - np.abs(c))
W = 5.0625  # root square [-W, W]^2 contains the disc
lv = 7; n = 2**lv
ix, iy = [v.ravel() for v in np.meshgrid(np.arange(n), np.arange(n), indexing="ij")]
t0 = time.time(); dmax = 0; tot = 0
while len(ix):
    size = 2*W/2**lv; rho = size/np.sqrt(2)
    c = (-W + (ix + 0.5)*size) + 1j*(-W + (iy + 0.5)*size)
    inside = np.abs(c) - rho < R_F
    d = np.where(inside, dist(c), -1.0)
    dmax = max(dmax, d.max()); tot += len(ix)
    keep = inside & (d + rho > THR)
    print("level %2d size %.3e cells %7d  remaining %6d  max d(centre) so far 1%+.3e" % (lv, size, len(ix), keep.sum(), dmax - 1), flush=True)
    if lv >= 46: break
    ix, iy = ix[keep], iy[keep]
    if len(ix) > 3_000_000: print("too many cells"); break
    ix = np.concatenate([2*ix, 2*ix + 1, 2*ix, 2*ix + 1]); iy = np.concatenate([2*iy, 2*iy, 2*iy + 1, 2*iy + 1]); lv += 1
if len(ix) == 0: print("WHOLE DISC COVERED: no point of the disc has distance-to-boundary > %.17g (double precision), %d cells, %.0f s" % (THR, tot, time.time() - t0))
else:
    print("NOT COVERED; remaining cell centres:"); 
    for z in c[keep][:40]: print("   ", z, abs(z))
