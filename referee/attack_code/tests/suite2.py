from mut import run
from suite1 import *
NA = [("    assert delta >= 0\n", "    assert delta.upper() >= 0\n")]
EL = [("a_half = math.sqrt(1 + b * b) * (1 + 1e-12)", "a_half = math.sqrt(1 + b * b) * (1 + 1e-2)"), ("r = rho_of(b) * (1 - 1e-12)", "r = rho_of(b) * (1 - 1e-2)")]
tests2 = [
 ("x_tq_1e-3", NA + S(TQ, "TOL_QUAD = 1e-3\n"), None), ("x_tq_1e-8", NA + S(TQ, "TOL_QUAD = 1e-8\n"), None), ("x_tq_1e2", NA + S(TQ, "TOL_QUAD = 1e2\n"), None),
 ("x_M600", NA + S(MLOOP, "        m_nodes = 600\n"), None),
 ("x_Nmin", NA + S(NLOOP, "            n_q = deg // 2 + 2\n"), None),
 ("x_tight2", NA + S(TI, "TOL_INTERP = 1e-30\n") + S(TQ, "TOL_QUAD = 1e-25\n"), None),
 ("x_el", NA + EL, None),
 ("x_prec53", NA + EL + S(PRECL, "PREC = 53\n"), None),
 ("x_prec40", NA + EL + S(PRECL, "PREC = 40\n"), None), ("x_prec30", NA + EL + S(PRECL, "PREC = 30\n"), None),
 ("x_prec24", NA + EL + S(PRECL, "PREC = 24\n"), None),
 ("x_prec20", NA + EL + S(PRECL, "PREC = 20\n"), None), ("x_prec12", NA + EL + S(PRECL, "PREC = 12\n"), None), ("x_prec6", NA + EL + S(PRECL, "PREC = 6\n"), None),
 ("x_prec40_tight", NA + EL + S(PRECL, "PREC = 40\n") + S(TI, "TOL_INTERP = 1e-14\n"), None),
 ("x_prec60_tight", NA + EL + S(PRECL, "PREC = 60\n") + S(TI, "TOL_INTERP = 1e-14\n"), None),
]
if __name__ == "__main__":
    only = sys.argv[2:]
    for d, n in (("ONE", 32), ("R6", 48)):
        if sys.argv[1] not in ("both", d): continue
        for tag, subs, de in tests2:
            if only and tag not in only: continue
            run(tag, subs, d, n, density_edit=de, timeout=240)
