import sys, time
import numpy as np
from bie import conformal_radius
from dom import R_CLAIM, claimed_slits

d2r = np.pi / 180
which = sys.argv[1]
ns = [int(v) for v in sys.argv[2].split(",")]
p = int(sys.argv[3]) if len(sys.argv) > 3 else 16

if which == "cal2":
    R = 3.0
    sl = [(k * 120 * d2r, 1.0, 0.0) for k in range(3)] + [((60 + k * 120) * d2r, 2.0, 0.0) for k in range(3)]
    ref = 1.5207430363125
elif which == "cal3":
    R = 5.99717
    sl = [(k * 120 * d2r, 1.0, 0.0) for k in range(3)] + [((60 + k * 120) * d2r, 1.99998, 0.0) for k in range(3)]
    for k in range(3):
        for sg in (1, -1):
            sl.append(((k * 120 + sg * 41.47240) * d2r, 3.87035, 0.0))
            sl.append(((k * 120 + sg * 19.61578) * d2r, 2.84221, 0.0))
    ref = 1.5235139885
elif which == "claim":
    R = R_CLAIM
    sl = [(al, tip, b) for (c, al, tip, b) in claimed_slits()]
    ref = 1.523673647065
elif which == "claim_nobend":      # diagnostic: same slits, bend ignored
    R = R_CLAIM
    sl = [(al, tip, 0.0) for (c, al, tip, b) in claimed_slits()]
    ref = 1.523673647065
elif which == "claim_flip":        # diagnostic: sign of bend flipped
    R = R_CLAIM
    sl = [(al, tip, -b) for (c, al, tip, b) in claimed_slits()]
    ref = 1.523673647065
for n in ns:
    t0 = time.time()
    r, _ = conformal_radius(R, sl, n=n, p=p)
    print(f"{which} n={n} p={p} N={len(sl)*n*p}: r={r:.13f}  r-ref={r-ref:+.3e}  ({time.time()-t0:.1f}s)", flush=True)
