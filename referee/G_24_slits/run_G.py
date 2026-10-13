import sys, time
import numpy as np
from bie import conformal_radius, make_edges
from dom import R_G, slits_G
CL = 1.523726248805
sl = [t[1:] for t in slits_G()]
cfgs = [dict(n=4, p=16), dict(n=6, p=16), dict(n=8, p=16), dict(n=12, p=16), dict(n=16, p=16),
        dict(n=8, p=10), dict(n=8, p=12), dict(n=8, p=20), dict(n=6, p=24),
        dict(n=10, p=16, m=2), dict(n=10, p=16, m=4), dict(n=10, p=16, adapt=1.5), dict(n=6, p=16, m=3, adapt=3.0),
        dict(n=10, p=16, power=4), dict(n=10, p=16, power=1),
        dict(p=16, edges=make_edges(10, 1.0, 2.0)), dict(p=16, edges=make_edges(10, 1.5, 1.0)),
        dict(p=16, edges=np.array([0, .03, .1, .2, .35, .5, .62, .75, .85, .93, .98, 1.0]))]
for c in cfgs:
    t0 = time.time()
    r, info = conformal_radius(R_G, sl, return_all=True, **c)
    d = {k: v for k, v in c.items() if k != "edges"}
    if "edges" in c: d["edges"] = "graded%d" % (len(c["edges"]) - 1)
    print(f"{d}: r={r:.14f} r-claim={r-CL:+.3e} min sig={info['sig'].min():+.2e} nadapt={info['nadapt']} ({time.time()-t0:.0f}s)", flush=True)
