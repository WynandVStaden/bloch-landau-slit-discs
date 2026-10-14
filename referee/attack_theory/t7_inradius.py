"""Independent floating-point look at the inradius (own code): sup over the sector of dist(p, K u circle)."""
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import minimize
from ubuild import *
h = 2e-5
pts = []
for j in range(NS):
    n = int((R-R*np.exp(-SL[j][1]))*1.7/h)+2
    pts.append(np.exp(wj(j, -np.linspace(0,1,n))))
pts.append(R*np.exp(1j*np.linspace(0,2*np.pi,int(2*np.pi*R/h))))
B = np.concatenate(pts); tree = cKDTree(np.column_stack([B.real,B.imag]))
print("boundary samples:", len(B), " max gap:", max(np.abs(np.diff(p)).max() for p in pts))
g = 0.004
X,Y = np.meshgrid(np.arange(0,R,g), np.arange(0,R,g)); P = (X+1j*Y).ravel()
P = P[(np.abs(P)<R)&(np.angle(P)<=np.pi/3+0.02)]
d,_ = tree.query(np.column_stack([P.real,P.imag]))
print("grid points:", len(P), " max grid dist:", d.max(), " => sup <= %.4f by Lipschitz" % (d.max()+g/np.sqrt(2)))
cand = P[d > 0.99]; dc = d[d>0.99]
seen = {}
for p,dd in zip(cand,dc):
    key = (round(p.real/0.1), round(p.imag/0.1))
    if key not in seen or seen[key][1] < dd: seen[key] = (p,dd)
f = lambda v: -tree.query([v[0],v[1]])[0]
res = []
for p,dd in seen.values():
    best = None
    r = minimize(f, [p.real,p.imag], method="Nelder-Mead", options={"xatol":1e-11,"fatol":1e-13,"initial_simplex":[[p.real,p.imag],[p.real+g,p.imag],[p.real,p.imag+g]], "maxiter":4000})
    res.append((-r.fun, r.x[0], r.x[1]))
res.sort(reverse=True)
uniq = []
for v in res:
    if all(abs(complex(v[1],v[2])-complex(u[1],u[2]))>1e-3 for u in uniq): uniq.append(v)
print("local maxima in the sector (sample-based distance, overestimates by <~ 1e-10):")
for v in uniq[:14]: print("   radius %.10f at (%.6f, %.6f)" % v)
print("largest:", uniq[0][0], " (claim: <= 1 + 1e-9)")
