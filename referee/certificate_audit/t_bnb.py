"""Does every ball built in bound_abs contain its float box? (exact rational check, G plan)"""
import json, sys
from fractions import Fraction
import numpy as np
import aud
cc = aud.cc
from flint import arb, acb
stats = {"boxes": 0, "viol": 0, "worst": Fraction(0), "minw": [9, 9, 9, 9], "calls": 0, "accepted": 0, "maxratio_hits": 0}
def bound_abs(func, box, scale, ratio=4.0, min_width=1e-7, max_boxes=400000):
    stats["calls"] += 1
    lo_all, hi_all = None, None
    stack, count = [box], 0
    while stack:
        bx = stack.pop(); count += 1
        if count > max_boxes: return None
        balls = [arb((lo + hi) / 2, (hi - lo) / 2 * (1 + 2.0**-30)) if hi > lo else arb(lo) for lo, hi in bx]
        for d, ((lo, hi), b) in enumerate(zip(bx, balls)):
            stats["boxes"] += 1
            m, r = Fraction(float(b.mid())), Fraction(float(b.rad()))
            assert arb(float(b.mid())) == b.mid() and (arb(float(b.rad())) == b.rad())
            gap = max(m - r - Fraction(lo), Fraction(hi) - (m + r))
            if gap > 0:
                stats["viol"] += 1; stats["worst"] = max(stats["worst"], gap)
            if hi > lo: stats["minw"][d] = min(stats["minw"][d], hi - lo)
        v = func(acb(balls[0], balls[1]), acb(balls[2], balls[3]))
        lo, hi = v.abs_lower(), v.abs_upper()
        if lo > 0 and hi < ratio * lo:
            stats["accepted"] += 1
            lo_all = lo if lo_all is None or lo < lo_all else lo_all
            hi_all = hi if hi_all is None or hi > hi_all else hi_all
            continue
        widths = [(h - l) * s for (l, h), s in zip(bx, scale)]
        d = int(np.argmax(widths))
        if widths[d] < min_width: return None
        mid = (bx[d][0] + bx[d][1]) / 2
        for part in ((bx[d][0], mid), (mid, bx[d][1])):
            child = list(bx); child[d] = part; stack.append(child)
    return lo_all, hi_all
cc.bound_abs = bound_abs
name = sys.argv[1] if len(sys.argv) > 1 else "G"
geo = cc.Geometry(name)
coef = {int(k): v for k, v in json.load(open(f"../../certify/conformal_certificate_{name}.json"))["coefficients"].items()}
pl = cc.plan(geo, coef)
print({k: (float(v) if isinstance(v, Fraction) else v) for k, v in stats.items()})
json.dump({str(i): {"b": p["b"], "rho": p["rho"], "B": p["B"], "M": p["M"], "interp_err": p["interp_err"], "quad": {str(j): q for j, q in p["quad"].items()}} for i, p in pl.items()}, open(f"plan_{name}.json", "w"))
