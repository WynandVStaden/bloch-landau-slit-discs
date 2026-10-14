# Part A: whole-disc search for local maxima of distance-to-boundary (double precision stage)
from geom import *
from scipy.spatial import cKDTree
from scipy import ndimage
from scipy.optimize import minimize
import pickle, time
t0 = time.time()
pts = []
for k, sl in enumerate(SL_F):
    pts.append(gam_f(sl, np.linspace(-sl['a'], 0, 40001)))
pts.append(R_F*np.exp(1j*np.linspace(0, 2*np.pi, 400001)))
P = np.concatenate(pts); tree = cKDTree(np.c_[P.real, P.imag])
h = 0.004; n = int(R_F/h) + 2
xs = np.arange(-n, n + 1)*h
D = np.zeros((2*n + 1, 2*n + 1), dtype=np.float32)
for i, y in enumerate(xs):
    q = np.c_[xs, np.full_like(xs, y)]
    d, _ = tree.query(q)
    d[np.hypot(xs, y) >= R_F] = 0
    D[i] = d
print("grid done", D.shape, time.time() - t0, "max grid d", D.max(), flush=True)
THR = 0.99
mask = D > THR
lab, nl = ndimage.label(mask, structure=np.ones((3, 3)))
print("components of {d>%.3f}: %d" % (THR, nl), flush=True)
def negd(v): return -dist_f(complex(v[0], v[1]))[0]
res = []
for L in range(1, nl + 1):
    ii, jj = np.nonzero(lab == L)
    vals = D[ii, jj]
    # starts: best grid point + grid local maxima + random points of component
    order = np.argsort(-vals)
    starts = [order[0]] + list(order[1:: max(1, len(order)//12)])
    best = None; found = []
    for o in starts:
        x0 = np.array([xs[jj[o]], xs[ii[o]]])
        for rep in range(3):
            r_ = minimize(negd, x0, method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-15, maxiter=4000, initial_simplex=x0 + (10.0**(-2 - 2*rep))*np.array([[0, 0], [1, 0], [0, 1]])))
            x0 = r_.x
        found.append((-r_.fun, r_.x[0], r_.x[1]))
    res.append(dict(label=L, size=len(ii), gridmax=float(vals.max()), cx=float(xs[jj].mean()), cy=float(xs[ii].mean()), found=found))
    f = sorted(found, reverse=True)
    print(L, "size", len(ii), "gridmax %.5f" % vals.max(), "opt: " + " ; ".join("%.12f @ (%.6f,%.6f)" % t for t in f[:3]), "spread", f[0][0] - f[-1][0], flush=True)
pickle.dump(res, open("search_res.pkl", "wb"))
print("time", time.time() - t0)
