import sys, io, contextlib
from fractions import Fraction
import numpy as np
sys.path.insert(0, "../../certify")
import inradius_certificate as ic
d = np.load("../../certify/inradius_leaves_G.npz")
cells, pars, eps = d["cells"], d["pars"], float(d["eps"])
print("stored leaves", cells.shape, "eps", repr(eps), "levels", cells[:, 0].min(), cells[:, 0].max())
inrange = ((cells[:, 1] >= 0) & (cells[:, 2] >= 0) & (cells[:, 1] < 2 ** cells[:, 0]) & (cells[:, 2] < 2 ** cells[:, 0])).all()
print("all cells inside root square:", bool(inrange), " codes:", np.bincount(cells[:, 3]), " kinds among covered:", np.bincount(cells[cells[:, 3] == 2, 4]))
big_r = Fraction(ic.DESIGNS["G"][0]); slits = ic.slit_table(ic.DESIGNS["G"][1])
ic.no_overlap(cells); print("stored G leaves verify:", ic.verify(big_r, slits, Fraction(eps), cells, pars))
# independent tiling check: rasterise coarse levels + recursive descent
leaf = {tuple(c[:3]) for c in cells.tolist()}
def covered(lv, x, y, depth=0):
    if (lv, x, y) in leaf: return True
    if lv >= 40: return False
    return all(covered(lv + 1, 2 * x + dx, 2 * y + dy) for dx in (0, 1) for dy in (0, 1))
print("independent descent: every level-6 cell is a leaf or fully subdivided into leaves:", all(covered(6, x, y) for x in range(64) for y in range(64)))
# attacks: the same leaves against perturbed domains must be rejected
def attack(tag, design, e=eps):
    sl = ic.slit_table(design[1]); br = Fraction(design[0])
    try:
        ic.verify(br, sl, Fraction(e), cells, pars); print(f"  {tag}: ACCEPTED")
    except AssertionError as ex:
        print(f"  {tag}: rejected ({ex})")
G = ic.DESIGNS["G"]
print("stored leaves vs perturbed claims:")
attack("eps = 0.9e-9", G, 0.9e-9)
attack("eps = 0", G, 0.0)
cl = list(G[1]); cl[0] = ("0", "1.00000001", "0", "0"); attack("class-0 tip 1 -> 1.00000001", (G[0], cl))
cl = list(G[1]); cl[4] = (cl[4][0], "4.3428742", cl[4][2], cl[4][3]); attack("class-4 tip moved 1.5e-8", (G[0], cl))
attack("R + 1e-8", ("5.001117638756668", G[1]))
# full pipeline on perturbed domains: must fail to certify 1e-9
for tag, des in [("G itself", G), ("class-0 tip 1.00000001", (G[0], [("0", "1.00000001", "0", "0")] + list(G[1][1:]))),
                 ("class-2 b1 + 1e-7", (G[0], list(G[1][:2]) + [(G[1][2][0], G[1][2][1], "0.9026922410641155", G[1][2][3])] + list(G[1][3:])))]:
    sl = ic.slit_table(des[1]); br = Fraction(des[0])
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            c2, p2 = ic.search(float(br), sl, 1e-9, max_level=36)
            ic.no_overlap(c2); res = ic.verify(br, sl, Fraction(1e-9), c2, p2)
        print(f"  full run, {tag}: CERTIFIED {res}")
    except (SystemExit, AssertionError) as ex:
        print(f"  full run, {tag}: not certified ({ex})")
