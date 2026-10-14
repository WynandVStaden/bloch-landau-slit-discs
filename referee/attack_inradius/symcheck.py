from geom import *
# symmetry: rotation by 2pi/3 and reflection z->conj(z) should permute slits
def key(sl): return (mp.nstr(mp.fmod(sl['th'] + 4*mp.pi, 2*mp.pi), 30), mp.nstr(sl['t'], 30), mp.nstr(sl['b1'], 30), mp.nstr(sl['b2'], 30))
def canon(th, t, b1, b2):
    th = mp.fmod(th + 8*mp.pi, 2*mp.pi)
    if abs(th - 2*mp.pi) < mpf(10)**-40: th = mpf(0)
    return (mp.nstr(th, 30), mp.nstr(t, 30), mp.nstr(b1, 30), mp.nstr(b2, 30))
S = set(canon(s['th'], s['t'], s['b1'], s['b2']) for s in SL_MP)
print("n slits", len(SL_MP), "distinct", len(S))
rot = set(canon(s['th'] + 2*mp.pi/3, s['t'], s['b1'], s['b2']) for s in SL_MP)
ref = set(canon(-s['th'], s['t'], -s['b1'], -s['b2']) for s in SL_MP)
print("rot invariant", rot == S, "reflection invariant", ref == S)
# pointwise check
import random
mx = 0
for k, sl in enumerate(SL_MP):
    for _ in range(5):
        s = -sl['a']*mpf(random.random()); z = gam(sl, s)
        for w in (z*mp.expj(2*mp.pi/3), mp.conj(z)):
            d = min(near_mp(j, w)[0] for j in range(24)); mx = max(mx, d)
print("max dist of image point to slit set", mp.nstr(mx, 5))
# min modulus of slit points, tips
for k, sl in enumerate(SL_MP):
    if sl['r'] == 0: print(k, sl['cls'], sl['eps'], "tip", mp.nstr(gam(sl, -sl['a']), 20), "a", mp.nstr(sl['a'], 20), "arg tip", mp.nstr(mp.arg(gam(sl, -sl['a'])), 15))
print("min tip modulus - 1:", mp.nstr(min(s['t'] for s in SL_MP) - 1, 5))
# pairwise min distances between slits (double)
import itertools
from scipy.spatial import cKDTree
for i, j in itertools.combinations(range(24), 2):
    if SL_F[i]['r'] != 0: continue
    A = np.c_[SAMP[i].real, SAMP[i].imag]; B = np.c_[SAMP[j].real, SAMP[j].imag]
    d, _ = cKDTree(B).query(A)
    if d.min() < 1.5: print("slits", i, (SL_F[i]['cls'], SL_F[i]['eps']), j, (SL_F[j]['cls'], SL_F[j]['r'], SL_F[j]['eps']), "min dist %.6f at |z|=%.4f" % (d.min(), abs(SAMP[i][d.argmin()])))
