"""The domains reported in the note; recomputes their ratios and draws them.

A slit is, in the coordinate w = log z,
    w(x) = log R + s + i (alpha + b1 s^2 + b2 s^4),   s = a x,  -1 <= x <= 0,  a = log(R / tip),
and each class (alpha, tip, b1, b2) stands for the slit and its images under
rotation by 120 degrees and complex conjugation.

Usage: python designs.py
"""
import numpy as np

import bloch_curved2 as b
from bloch_solver import conformal_radius

PROVED_2009 = 0.6563937  # Carroll and Ortega-Cerda, explicit formula
NUMERICAL_2010 = 0.656319277272  # Ortega-Cerda and Pridhnani, no rigorous error analysis

DESIGNS = {
    "A: 18 straight slits": {
        "R": 5.23368,
        "classes": [(0.0, 1.0, 0.0, 0.0), (np.pi / 3, 2.0, 0.0, 0.0),
                    (0.80248, 3.74887, 0.0, 0.0), (0.32992, 2.83043, 0.0, 0.0)],
    },
    "B: 18 slits, 12 curved": {
        "R": 4.747749358688023,
        "classes": [(0.0, 1.0, 0.0, 0.0), (np.pi / 3, 1.9999999999999036, 0.0, 0.0),
                    (0.12393306687411461, 2.8051388295325297, 0.6738078607105484, 0.0),
                    (0.7818321423345288, 3.710751140554841, -0.6978713930533904, 0.0)],
    },
    "G: 24 slits, 18 curved": {
        "R": 5.001117628756668,
        "classes": [(0.0, 1.0, 0.0, 0.0), (np.pi / 3, 1.9999999974092333, 0.0, 0.0),
                    (0.029930188727091613, 2.789342409403121, 0.9026921410641155, -0.3301459494660587),
                    (1.0157094972537686, 3.473381248030611, -1.2656124700253997, 0.9794269477754911),
                    (0.5757776442359243, 4.342874184977691, 0.03225521024585182, 2.9360081049550737)],
    },
}


def evaluate(design, n=96):
    big_r, classes = design["R"], design["classes"]
    inr = b.inradius(big_r, classes, grid=(400, 201), per_slit=800)
    conf = conformal_radius(big_r, b.slit_objects(big_r, classes), n)
    return inr, conf, inr / conf


def draw(design, ax, title):
    import matplotlib.patches as patches
    big_r, classes = design["R"], design["classes"]
    ang, tip, b1, b2 = b.expand(classes)
    a = np.log(big_r / tip)
    t = np.linspace(0, 2 * np.pi, 720)
    ax.plot(big_r * np.cos(t), big_r * np.sin(t), color="black", lw=1.1)
    x = -np.linspace(0, 1, 400)
    for k in range(len(ang)):
        z = b.curve(x, big_r, ang[k], a[k], b1[k], b2[k])
        ax.plot(z.real, z.imag, color="black", lw=1.1)
    ax.add_patch(patches.Circle((0, 0), 1.0, fill=False, color="#c44e52", lw=1.0, ls="--"))
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.axis("off")


if __name__ == "__main__":
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4))
    for (name, design), ax in zip(DESIGNS.items(), axes):
        inr, conf, ratio = evaluate(design)
        print(f"{name}: R = {design['R']}, inradius {inr:.10f}, conformal radius {conf:.12f}, ratio {ratio:.9f}; "
              f"below 0.6563937 by {PROVED_2009 - ratio:+.3e}, below 0.6563193 by {NUMERICAL_2010 - ratio:+.3e}")
        draw(design, ax, f"{name}\nratio {ratio:.7f}")
    fig.suptitle("Slit discs with inradius 1 (dashed: a largest inscribed disc)", fontsize=11)
    fig.tight_layout()
    fig.savefig("domains.png", dpi=150)
    print("figure written to domains.png")
