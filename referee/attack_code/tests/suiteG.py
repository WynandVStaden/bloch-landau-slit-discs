from mut import run
from suite1 import *
from suite2 import NA, EL
def c1cls(cls, f):
    def e(d): d["coefficients"][cls][1] *= f
    return e
def setc(cls, n, f):
    def e(d): d["coefficients"][cls][n] *= f
    return e
testsG = [
 ("G_c1dn_cls0", [], c1cls("0", 1 - 1e-6)), ("G_c1up_cls0", [], c1cls("0", 1 + 1e-6)),
 ("G_c1dn_cls4", [], c1cls("4", 1 - 1e-4)),
 ("G_c63_cls2", [], setc("2", 63, 100.0)),
 ("G_saf09", S(SAF, "SAFETY = 0.9 "), None), ("G_saf099", S(SAF, "SAFETY = 0.99 "), None),
 ("G_BUG_drop5", S(QLOOP, QLOOP + "            if j == 5: continue\n"), None),
 ("G_BUG_drop23", S(QLOOP, QLOOP + "            if j == 23: continue\n"), None),
 ("G_BUG_N8", S(SRCN, SRCN.replace("n_q)", "max(1, n_q // 8))")), None),
 ("G_BUG_Nhalf", S(SRCN, SRCN.replace("n_q)", "max(1, n_q // 2))")), None),
 ("G_BUG_logflip", S(LOGP, LOGP.replace("f +=", "f -=")), None),
 ("G_prec53", S(PRECL, "PREC = 53\n"), None),
 ("G_prec40el", NA + EL + S(PRECL, "PREC = 40\n"), None),
 ("G_ti1e-11", S(TI, "TOL_INTERP = 1e-11\n"), None),
 ("G_ratio1e3", S("ratio=4.0", "ratio=1e3"), None),
]
if __name__ == "__main__":
    only = sys.argv[1:]
    for tag, subs, de in testsG:
        if only and tag not in only: continue
        run(tag, subs, "G", 64, density_edit=de, timeout=1500)
