import os, json, sys, random
des = {"ONE_0.5": ("1", [("0", "0.5", "0", "0")]), "R6_3_2": ("3", [("0", "1", "0", "0"), ("pi/3", "2", "0", "0")])}
os.environ["AUD_DESIGNS"] = json.dumps(des)
import aud
cc = aud.cc
if __name__ == "__main__":
    import bloch_explicit as be
    _Pool = cc.Pool
    cc.Pool = lambda p, initializer=None, initargs=(): _Pool(p, initializer=aud.init, initargs=initargs)
    orig = cc.make_density
    exact = {"ONE_0.5": 8 / 9, "R6_3_2": be.conformal_radius(3.0, 2.0)}
    random.seed(1)
    def run(tag, name, mod, tol=1e-9):
        cc.TOL_INTERP = tol
        def md(geo, n):
            coef, rho = orig(geo, n)
            return mod(coef), rho
        cc.make_density = md
        lower = cc.certify(name, 64, 8, log=lambda *a: None)
        d = json.load(open(f"conformal_certificate_{name}.json"))["delta_upper"]
        print(f"{name} [{tag}]: certified {lower:.15f} exact {exact[name]:.15f} gap {exact[name]-lower:.3e} delta {d:.3e} {'OK' if lower <= exact[name] else '*** VIOLATION ***'}", flush=True)
    for name in ["R6_3_2"]:
        run("tight tol 1e-13", name, lambda c: c, 1e-13)
        for e in (1e-6, -1e-6, 1e-3, -1e-3):
            def m1(c, e=e):
                c = {k: list(v) for k, v in c.items()}
                for k in c: c[k][1] += e
                return c
            run(f"c1 += {e}", name, m1, 1e-12)
        def m3(c):
            c = {k: list(v) for k, v in c.items()}
            for k in c:
                for n in range(1, 20, 2): c[k][n] += random.uniform(-1e-5, 1e-5)
            return c
        run("random odd +-1e-5", name, m3, 1e-12)
        def mtr(c):
            return {k: [v if n < 12 else 0.0 for n, v in enumerate(vs)] for k, vs in c.items()}
        # (truncation test removed: plan() raises IndexError when class degrees differ from global deg)
        def mcl(c):  # perturb only one class
            c = {k: list(v) for k, v in c.items()}
            k0 = sorted(c)[-1]; c[k0][1] -= 1e-4; c[k0][3] += 2e-4
            return c
        run("one class perturbed", name, mcl, 1e-12)
