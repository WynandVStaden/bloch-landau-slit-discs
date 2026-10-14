# dist() extracted from quadcover.py (referee)
import numpy as np
from geom import *
from scipy.spatial import cKDTree
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
