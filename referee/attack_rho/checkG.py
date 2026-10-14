"""Interior sign/consistency check of the Green function: G>0 in D, and N-independence away from slits."""
import numpy as np
from geom import full_problem
from methodA import solve, eval_G
from methodC_wos import dist, full
R, classes = full_problem()
rng = np.random.default_rng(1)
z = R * np.sqrt(rng.random(4000)) * np.exp(2j * np.pi * rng.random(4000))
d = np.minimum(dist(z, full[0], full[1])[0], R - abs(z))
G = {}
for N in (192, 768):
    f, C, A, rhs = solve(R, classes, N)
    G[N] = eval_G(R, classes, N, f, C, z)
for N in G:
    neg = G[N] < 0
    print(f'N={N}: #negative={neg.sum()}, max dist-to-boundary among negative points={d[neg].max() if neg.any() else 0:.2e}, min G={G[N].min():.2e}')
far = d > 0.05
print('points with dist>0.05:', far.sum(), ' max |G_192-G_768| there:', np.abs(G[192] - G[768])[far].max(), ' min G there:', G[768][far].min())
