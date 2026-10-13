"""Rigorous upper bound for the inradius of a slit disc.

Claim proved by a successful run:  every disc contained in D has radius <= 1 + EPS.

D is the disc |z| < R minus the slits; it is invariant under rotation by 120
degrees and under conjugation, so it is enough to bound dist(p, boundary) for p
in the closed sector 0 <= arg p <= 60 degrees, |p| <= R.

Stage 1 (floating point, not trusted): a quadtree over the square [0, 8]^2.
A cell becomes a leaf when
   (a) it lies outside the disc, or
   (b) it lies outside the sector, or
   (c) a boundary point q ("witness") is found with all four corners of the
       cell within 1 + EPS of q. Then every p in the cell has
       dist(p, boundary) <= |p - q| <= 1 + EPS, since |p - q| is convex in p.
The leaves tile the square by construction (each split cell is replaced by its
four children).

Stage 2 (trusted): every leaf is re-checked.
   * Cell corners are dyadic rationals, held exactly as integers.
   * R and the slit parameters are the decimal numbers in the table, held
     exactly; the angle pi/3 is exact.
   * A witness is a point of a slit given by its parameter s (a double, hence
     an exact rational), a slit tip, or a point of the circle. Its coordinates
     are enclosed by ball arithmetic (Arb, via python-flint) and rounded to a
     dyadic rational q' with |q - q'| < 2^-69, and -a < s <= 0 is checked.
   * (a) is checked as xmin^2 + ymin^2 >= R^2, (b) as ymin^2 > 3 xmax^2,
     (c) as |corner - q'|^2 <= (1 + EPS - 2^-69)^2, all in exact integer
     arithmetic.

Usage: python inradius_certificate.py [design letter, default G] [EPS, default 1e-7]
"""
import sys
import time
from fractions import Fraction

import numpy as np
from flint import arb, ctx
from scipy.spatial import cKDTree

PI3 = "pi/3"
DESIGNS = {
    "B": ("4.747749358688023", [("0", "1", "0", "0"), (PI3, "1.9999999999999036", "0", "0"),
                                ("0.12393306687411461", "2.8051388295325297", "0.6738078607105484", "0"),
                                ("0.7818321423345288", "3.710751140554841", "-0.6978713930533904", "0")]),
    "G": ("5.001117628756668", [("0", "1", "0", "0"), (PI3, "1.9999999974092333", "0", "0"),
                                ("0.029930188727091613", "2.789342409403121", "0.9026921410641155", "-0.3301459494660587"),
                                ("1.0157094972537686", "3.473381248030611", "-1.2656124700253997", "0.9794269477754911"),
                                ("0.5757776442359243", "4.342874184977691", "0.03225521024585182", "2.9360081049550737")]),
}
ROOT = 8  # the quadtree covers [0, ROOT]^2
SCALE_BITS = 80
ETA_BITS = 69


def slit_table(classes):
    """All slits as (angle numerator data, tip, b1, b2) with exact entries: angle = 2 pi r / 3 + sgn * alpha."""
    out = []
    for alpha, tip, b1, b2 in classes:
        axial = alpha in ("0", PI3)
        for r in range(3):
            for sgn in ((1,) if axial else (1, -1)):
                out.append((r, sgn, alpha, Fraction(tip), sgn * Fraction(b1), sgn * Fraction(b2)))
    return out


def arb_frac(q):
    return arb(q.numerator) / arb(q.denominator)


def arb_angle(r, sgn, alpha):
    al = arb.pi() / 3 if alpha == PI3 else arb_frac(Fraction(alpha))
    return 2 * arb.pi() * r / 3 + sgn * al


def float_angle(r, sgn, alpha):
    return 2 * np.pi * r / 3 + sgn * (np.pi / 3 if alpha == PI3 else float(alpha))


# ---------------------------------------------------------------- stage 1

def build_samples(big_r, slits, spacing):
    pts, kind, idx, par = [], [], [], []
    for j, (r, sgn, alpha, tip, b1, b2) in enumerate(slits):
        a = np.log(big_r / float(tip))
        ang = float_angle(r, sgn, alpha)
        n = int(np.ceil((big_r - float(tip)) * 1.6 / spacing)) + 2
        s = -a * np.linspace(0, 1, n)[:-1] * (1 - 1e-12)  # strictly inside (-a, 0]
        z = big_r * np.exp(s + 1j * (ang + float(b1) * s**2 + float(b2) * s**4))
        pts.append(z), kind.append(np.zeros(len(s), int)), idx.append(np.full(len(s), j)), par.append(s)
        s_tip = -a
        pts.append(np.array([big_r * np.exp(s_tip + 1j * (ang + float(b1) * s_tip**2 + float(b2) * s_tip**4))]))
        kind.append(np.ones(1, int)), idx.append(np.full(1, j)), par.append(np.zeros(1))
    th = np.linspace(-0.3, np.pi / 3 + 0.3, int(np.ceil((np.pi / 3 + 0.6) * big_r / spacing)))
    pts.append(big_r * np.exp(1j * th)), kind.append(np.full(len(th), 2)), idx.append(np.zeros(len(th), int)), par.append(th)
    return np.concatenate(pts), np.concatenate(kind), np.concatenate(idx), np.concatenate(par)


def search(big_r, slits, eps, spacing=2e-4, max_level=44):
    pts, kind, idx, par = build_samples(big_r, slits, spacing)
    tree = cKDTree(np.column_stack([pts.real, pts.imag]))
    level = 6
    n = 2**level
    ix, iy = [v.ravel() for v in np.meshgrid(np.arange(n), np.arange(n), indexing="ij")]
    leaves = []  # (level, ix, iy, code, witness kind, witness slit, witness parameter)
    while len(ix):
        size = ROOT / 2**level
        x0, y0 = ix * size, iy * size
        outside_disc = np.hypot(x0, y0) > big_r + 1e-9
        outside_sector = y0 > np.sqrt(3) * (x0 + size) + 1e-9
        d0, k = tree.query(np.column_stack([x0 + size / 2, y0 + size / 2]))
        covered = d0 + size / np.sqrt(2) <= 1 + eps - 1e-10
        code = np.where(outside_disc, 0, np.where(outside_sector, 1, np.where(covered, 2, -1)))
        done = code >= 0
        for c in (0, 1, 2):
            m = code == c
            if m.any():
                leaves.append(np.column_stack([np.full(m.sum(), level), ix[m], iy[m], np.full(m.sum(), c),
                                               kind[k[m]], idx[k[m]]]).astype(np.int64))
                leaves.append(par[k[m]])
        ix, iy = ix[~done], iy[~done]
        print(f"  level {level:2d}: cell size {size:.2e}, {done.sum():8d} leaves, {len(ix):7d} cells to split", flush=True)
        if level == max_level and len(ix):
            raise SystemExit("search failed: cells remain at the maximum depth; EPS is too small for this domain")
        ix = np.concatenate([2 * ix, 2 * ix + 1, 2 * ix, 2 * ix + 1])
        iy = np.concatenate([2 * iy, 2 * iy, 2 * iy + 1, 2 * iy + 1])
        level += 1
    return np.concatenate(leaves[0::2]), np.concatenate(leaves[1::2])


# ---------------------------------------------------------------- stage 2

def to_fixed(x):
    """Integer n with |x - n / 2^SCALE_BITS| < 2^-(ETA_BITS + 1), for an Arb ball x."""
    assert x.rad() < arb(2) ** (-(ETA_BITS + 3))
    n = int((x.mid() * arb(2) ** SCALE_BITS).floor().unique_fmpz())
    return n


def witness_point(big_r, slits, kind, j, par):
    """Enclosure of a boundary point; raises if the parameter is not certified to lie on the slit."""
    r_ball = arb_frac(big_r)
    if kind == 2:
        th = arb_frac(Fraction(float(par)))
        return r_ball * th.cos(), r_ball * th.sin()
    r, sgn, alpha, tip, b1, b2 = slits[j]
    a = (r_ball / arb_frac(tip)).log()
    if kind == 1:
        s = -a
    else:
        s = arb_frac(Fraction(float(par)))
        assert par <= 0 and (s + a > 0) is True, "parameter not certified to lie in (-a, 0]"
    th = arb_angle(r, sgn, alpha) + arb_frac(b1) * s**2 + arb_frac(b2) * s**4
    rad = r_ball * s.exp()
    return rad * th.cos(), rad * th.sin()


def verify(big_r, slits, eps, cells, pars):
    ctx.prec = 160
    one = 1 << SCALE_BITS
    limit = one + int(eps * one) - (1 << (SCALE_BITS - ETA_BITS))  # (1 + eps - 2^-69), scaled, rounded down
    limit2 = limit * limit
    r2_num, r2_den = (big_r * big_r).numerator, (big_r * big_r).denominator
    cache, counts = {}, [0, 0, 0]
    # the leaves must tile the root square: total area check in exact arithmetic
    area_by_level = {int(lv): int((cells[:, 0] == lv).sum()) for lv in np.unique(cells[:, 0])}
    area = sum(Fraction(cnt, 4**lv) for lv, cnt in area_by_level.items())
    assert area == 1, f"leaves do not have the area of the root square (got {area})"
    for (lv, ix, iy, code, kind, j), par in zip(cells.tolist(), pars.tolist()):
        shift = SCALE_BITS + 3 - lv  # cell size = 8 / 2^lv = 2^(3 - lv)
        assert shift >= 0
        x0, y0 = ix << shift, iy << shift
        x1, y1 = x0 + (1 << shift), y0 + (1 << shift)
        if code == 0:
            assert (x0 * x0 + y0 * y0) * r2_den >= r2_num * one * one, "cell not outside the disc"
        elif code == 1:
            assert y0 * y0 > 3 * x1 * x1, "cell not outside the sector"
        else:
            key = (kind, j, par)
            if key not in cache:
                qx, qy = witness_point(big_r, slits, kind, j, par)
                cache[key] = (to_fixed(qx), to_fixed(qy))
            qx, qy = cache[key]
            for cx in (x0, x1):
                for cy in (y0, y1):
                    assert (cx - qx) ** 2 + (cy - qy) ** 2 <= limit2, "corner too far from the witness"
        counts[code] += 1
    return counts, len(cache)


def no_overlap(cells):
    """Leaves come from a quadtree, so equal area plus pairwise distinctness of ancestors gives a tiling; check no leaf is an ancestor of another."""
    seen = set(map(tuple, cells[:, :3].tolist()))
    assert len(seen) == len(cells)
    for lv, ix, iy in seen:
        l2, x, y = lv, ix, iy
        while l2 > 6:
            l2, x, y = l2 - 1, x >> 1, y >> 1
            assert (l2, x, y) not in seen, "a leaf contains another leaf"


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "G"
    eps = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-7
    big_r = Fraction(DESIGNS[name][0])
    slits = slit_table(DESIGNS[name][1])
    assert big_r < ROOT
    print(f"design {name}: R = {DESIGNS[name][0]}, {len(slits)} slits, EPS = {eps}")
    t = time.time()
    print("stage 1: quadtree search (floating point)")
    cells, pars = search(float(big_r), slits, eps)
    print(f"  {len(cells)} leaves in {time.time() - t:.0f} s")
    np.savez_compressed(f"inradius_leaves_{name}.npz", cells=cells, pars=pars, eps=eps)
    t = time.time()
    print("stage 2: exact verification of every leaf")
    no_overlap(cells)
    counts, n_wit = verify(big_r, slits, Fraction(eps), cells, pars)
    print(f"  outside disc {counts[0]}, outside sector {counts[1]}, covered by a witness {counts[2]} "
          f"({n_wit} distinct witnesses); {time.time() - t:.0f} s")
    print(f"CERTIFIED: inradius of design {name} <= 1 + {eps}")
