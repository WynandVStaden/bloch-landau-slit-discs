"""Domain definition, written from the referee brief only."""
import numpy as np

R_CLAIM = 4.747749358688023


def claimed_slits():
    s = []
    for k in range(3):
        s.append(("A", 2 * np.pi * k / 3, 1.0, 0.0))
    for k in range(3):
        s.append(("B", np.pi / 3 + 2 * np.pi * k / 3, 1.9999999999999036, 0.0))
    for k in range(3):
        for sg in (1, -1):
            s.append(("C", 2 * np.pi * k / 3 + sg * 0.12393306687411461,
                      2.8051388295325297, sg * 0.6738078607105484))
    for k in range(3):
        for sg in (1, -1):
            s.append(("D", 2 * np.pi * k / 3 + sg * 0.7818321423345288,
                      3.710751140554841, sg * (-0.6978713930533904)))
    return s


def zslit(x, R, alpha, tip, bend):
    """Point of slit at parameter x in [-1,0]."""
    a = np.log(R / tip)
    u = a * x
    return R * np.exp(u + 1j * (alpha + bend * u * u))


def dzslit(x, R, alpha, tip, bend):
    a = np.log(R / tip)
    u = a * x
    return zslit(x, R, alpha, tip, bend) * (a + 1j * 2 * bend * u * a)
