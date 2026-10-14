# second part of mutate_full.py; search depth capped at 34 (true domain needs 32) to keep refusals fast
# Full-pipeline mutation tests: run the authors' search + verify (from the COPY) on false claims.
import sys, io, contextlib, numpy as np
from fractions import Fraction
sys.path.insert(0, "run_copy")
import inradius_certificate_copy as P

def run(name, R, classes, eps_text, expect):
    big_r = Fraction(R); slits = P.slit_table(classes); buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            cells, pars = P.search(float(big_r), slits, float(eps_text), max_level=34)
            P.no_overlap(cells)
            counts, n_wit = P.verify(big_r, slits, Fraction(eps_text), cells, pars)
        res = "CERTIFIED (%d leaves, max level %d)" % (len(cells), cells[:, 0].max())
    except SystemExit as e:
        res = "refused: " + str(e)[:60]
    except BaseException as e:
        res = "refused: %s %s" % (type(e).__name__, str(e)[:60])
    bad = (res.startswith("CERTIFIED") != (expect == "CERTIFIED")) and expect != "?"
    print("%-58s eps=%-8s expect %-9s -> %s%s" % (name, eps_text, expect, res, "   ###### UNEXPECTED" if bad else ""), flush=True)

R0, C0 = P.DESIGNS["G"]; C0 = [tuple(c) for c in C0]
def rep(i, k, v):
    c = [list(x) for x in C0]; c[i][k] = v; return [tuple(x) for x in c]
run("R -> 5.0012", "5.0012", C0, "1e-9", "?")
run("class2 b1 -> 0.9027 (tip rotates by 2.7e-6 rad)", R0, rep(2, 2, "0.9027"), "1e-9", "?")
run("class2 b1 -> 0.9026 (tip rotates the other way)", R0, rep(2, 2, "0.9026"), "1e-9", "?")
run("class3 theta -> 1.0158", R0, rep(3, 0, "1.0158"), "1e-9", "?")
run("class3 theta -> 1.0156", R0, rep(3, 0, "1.0156"), "1e-9", "?")
run("class4 removed", R0, C0[:4], "1e-9", "refused")
run("class3 removed", R0, C0[:3] + C0[4:], "1e-9", "refused")
run("class2 removed", R0, C0[:2] + C0[3:], "1e-9", "refused")
run("class1 removed", R0, C0[:1] + C0[2:], "1e-9", "refused")
run("R -> 5.1", "5.1", C0, "1e-9", "?")
