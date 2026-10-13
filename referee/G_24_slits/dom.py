"""Domain definitions written from the referee brief only (slits with b1*s^2 + b2*s^4)."""
import numpy as np

R_G = 5.001117628756668
TWO3 = 2 * np.pi / 3


def slits_G():
    """list of (cls, alpha, tip, b1, b2)"""
    s = []
    for k in range(3):
        s.append(("A", TWO3 * k, 1.0, 0.0, 0.0))
    for k in range(3):
        s.append(("B", np.pi / 3 + TWO3 * k, 1.9999999974092333, 0.0, 0.0))
    for k in range(3):
        for sg in (1, -1):
            s.append(("C", TWO3 * k + sg * 0.029930188727091613, 2.789342409403121,
                      sg * 0.9026921410641155, sg * (-0.3301459494660587)))
    for k in range(3):
        for sg in (1, -1):
            s.append(("D", TWO3 * k + sg * 1.0157094972537686, 3.473381248030611,
                      sg * (-1.2656124700253997), sg * 0.9794269477754911))
    for k in range(3):
        for sg in (1, -1):
            s.append(("E", TWO3 * k + sg * 0.5757776442359243, 4.342874184977691,
                      sg * 0.03225521024585182, sg * 2.9360081049550737))
    return s


R_18 = 4.747749358688023
REF_18 = 1.5236736470653


def slits_18():
    s = []
    for k in range(3):
        s.append(("A", TWO3 * k, 1.0, 0.0, 0.0))
    for k in range(3):
        s.append(("B", np.pi / 3 + TWO3 * k, 1.9999999999999036, 0.0, 0.0))
    for k in range(3):
        for sg in (1, -1):
            s.append(("C", TWO3 * k + sg * 0.12393306687411461, 2.8051388295325297,
                      sg * 0.6738078607105484, 0.0))
    for k in range(3):
        for sg in (1, -1):
            s.append(("D", TWO3 * k + sg * 0.7818321423345288, 3.710751140554841,
                      sg * (-0.6978713930533904), 0.0))
    return s


def theta(s, alpha, b1, b2):
    return alpha + b1 * s * s + b2 * s ** 4


def zslit(x, R, alpha, tip, b1, b2):
    """x in [-1,0]; x=-1 tip, x=0 on circle."""
    a = np.log(R / tip)
    s = a * x
    return R * np.exp(s + 1j * theta(s, alpha, b1, b2))


def dzslit(x, R, alpha, tip, b1, b2):
    a = np.log(R / tip)
    s = a * x
    return zslit(x, R, alpha, tip, b1, b2) * a * (1 + 1j * (2 * b1 * s + 4 * b2 * s ** 3))
