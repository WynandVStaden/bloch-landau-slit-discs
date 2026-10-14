from mut import run
import sys
S = lambda a, b: [(a, b)]
SAF = "SAFETY = 0.5 "
PRECL = "PREC = 100\n"
TI = "TOL_INTERP = 1e-9\n"
TQ = "TOL_QUAD = 1e-13\n"
MLOOP = "        while not interp_remainder(m_nodes) < TOL_INTERP:\n            m_nodes += 2\n"
NLOOP = "            while not quad_remainder(n_q) < TOL_QUAD:\n                n_q += 1\n"
QLOOP = "        for j, (n_q, err) in quad.items():\n"
LOGP = "                f += arb(c) * (arb.pi() / n) * arb.cos_pi_fmpq(fmpq(n * (2 * k + 1), 2 * m_nodes))"
SRCN = "            ws, ps, ys = _source_nodes(j, n_q)\n"
def scale_c1(f):
    def e(d):
        for c in d["coefficients"]: d["coefficients"][c][1] *= f
    return e
def scale_all(f):
    def e(d):
        for c in d["coefficients"]: d["coefficients"][c] = [v * f for v in d["coefficients"][c]]
    return e
def even(d):
    d["coefficients"]["0"][2] = 1e-300
def trunc(nkeep):
    def e(d):
        for c in d["coefficients"]:
            d["coefficients"][c] = [v if n < nkeep else 0.0 for n, v in enumerate(d["coefficients"][c])]
    return e
def nan(d): d["coefficients"]["0"][3] = float("nan")
tests = [
 ("base", [], None),
 ("tight", S(TI, "TOL_INTERP = 1e-14\n") + S(TQ, "TOL_QUAD = 1e-17\n"), None),
 ("c1_up", [], scale_c1(1.001)), ("c1_dn", [], scale_c1(0.999)), ("c1_dn_tight", S(TI, "TOL_INTERP = 1e-14\n"), scale_c1(1 - 1e-9)),
 ("c1_up_tight", S(TI, "TOL_INTERP = 1e-14\n"), scale_c1(1 + 1e-9)),
 ("all_09", [], scale_all(0.9)), ("all_11", [], scale_all(1.1)), ("all_neg", [], scale_all(-1.0)), ("trunc8", [], trunc(8)), ("trunc2", [], trunc(2)),
 ("even", [], even), ("nan", [], nan),
 ("saf09", S(SAF, "SAFETY = 0.9 "), None), ("saf099", S(SAF, "SAFETY = 0.99 "), None), ("saf3", S(SAF, "SAFETY = 3.0 "), None),
 ("ti_1e-2", S(TI, "TOL_INTERP = 1e-2\n"), None), ("ti_1e-5", S(TI, "TOL_INTERP = 1e-5\n"), None),
 ("tq_1e-3", S(TQ, "TOL_QUAD = 1e-3\n"), None), ("tq_1e-6", S(TQ, "TOL_QUAD = 1e-6\n"), None), ("tq_1e2", S(TQ, "TOL_QUAD = 1e2\n"), None),
 ("prec53", S(PRECL, "PREC = 53\n"), None), ("prec40", S(PRECL, "PREC = 40\n"), None), ("prec30", S(PRECL, "PREC = 30\n"), None), ("prec20", S(PRECL, "PREC = 20\n"), None), ("prec10", S(PRECL, "PREC = 10\n"), None), ("prec4", S(PRECL, "PREC = 4\n"), None),
 ("ratio1e6", S("ratio=4.0", "ratio=1e6"), None), ("ratio1.2", S("ratio=4.0", "ratio=1.2"), None),
 ("M4", S(MLOOP, "        m_nodes = 4\n"), None), ("M2", S(MLOOP, "        m_nodes = 2\n"), None), ("M20", S(MLOOP, "        m_nodes = 20\n"), None), ("M600", S(MLOOP, "        m_nodes = 600\n"), None),
 ("Nmin", S(NLOOP, "            n_q = deg // 2 + 2\n"), None), ("N500", S(NLOOP, "            n_q = 500\n"), None),
 # injected bugs
 ("BUG_drop_self", S(QLOOP, QLOOP + "            if j == 0: continue\n"), None),
 ("BUG_drop_last", S(QLOOP, QLOOP + "            if j == max(quad): continue\n"), None),
 ("BUG_logflip", S(LOGP, LOGP.replace("f +=", "f -=")), None),
 ("BUG_Nhalf", S(SRCN, SRCN.replace("n_q)", "max(1, n_q // 2))")), None),
 ("BUG_N8", S(SRCN, SRCN.replace("n_q)", "max(1, n_q // 8))")), None),
 ("BUG_N16", S(SRCN, SRCN.replace("n_q)", "max(1, n_q // 16))")), None),
 ("BUG_noax", S("        f = -geo.a[i] * x\n", "        f = -geo.a[i] * x * (1 - 1e-7)\n"), None),
 ("BUG_kernel", S("                    acc -= abs(geo.self_kernel(i, x, y)).log() * p\n", "                    acc -= abs(geo.self_kernel(i, x, y) * (1 + 1e-7)).log() * p\n"), None),
]
if __name__ == '__main__':
    only = sys.argv[2:] 
    for d, n in (("ONE", 32), ("R6", 48)):
        if sys.argv[1] not in ("both", d): continue
        for tag, subs, de in tests:
            if only and tag not in only: continue
            run(tag, subs, d, n, density_edit=de, timeout=300)
