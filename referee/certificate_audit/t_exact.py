import os, json, sys
des = {}
for s in ["0.1", "0.5", "0.9", "0.99"]:
    des["ONE_" + s] = ("1", [("0", s, "0", "0")])
for R, t in [("3", "2"), ("5", "1.5"), ("2.5", "2.4")]:
    des[f"R6_{R}_{t}"] = (R, [("0", "1", "0", "0"), ("pi/3", t, "0", "0")])
os.environ["AUD_DESIGNS"] = json.dumps(des)
import aud
cc = aud.cc
import multiprocessing as mp
if __name__ == "__main__":
    import bloch_explicit as be
    # make workers use aud.init so that they see the extra designs
    import functools
    _Pool = cc.Pool
    cc.Pool = lambda p, initializer=None, initargs=(): _Pool(p, initializer=aud.init, initargs=initargs)
    for name in des:
        try:
            lower = cc.certify(name, 64, 8, log=lambda *a: None)
        except Exception as e:
            print(name, "FAILED TO RUN", repr(e)); continue
        if name.startswith("ONE"):
            s = float(name.split("_")[1]); exact = 4 * s / (1 + s) ** 2
        else:
            _, R, t = name.split("_"); exact = be.conformal_radius(float(R), float(t))
        print(f"{name}: certified {lower:.15f} exact {exact:.15f} gap {exact-lower:.3e} {'OK' if lower <= exact else '*** VIOLATION ***'}", flush=True)
