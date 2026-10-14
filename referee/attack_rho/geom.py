"""Geometry of the slit disc (independent re-implementation from the problem statement)."""
import numpy as np

R0 = 5.001117628756668
CLASSES0 = [
    (0.0, 1.0, 0.0, 0.0),
    (np.pi / 3, 1.9999999974092333, 0.0, 0.0),
    (0.029930188727091613, 2.789342409403121, 0.9026921410641155, -0.3301459494660587),
    (1.0157094972537686, 3.473381248030611, -1.2656124700253997, 0.9794269477754911),
    (0.5757776442359243, 4.342874184977691, 0.03225521024585182, 2.9360081049550737),
]


def full_problem():
    """Return (R, classes) with classes = list of dict(theta,t,b1,b2, images=[(rot, conj)])."""
    cl = []
    for i, (th, t, b1, b2) in enumerate(CLASSES0):
        if i < 2:
            imgs = [(2 * np.pi * r / 3, False) for r in range(3)]
        else:
            imgs = [(2 * np.pi * r / 3, c) for r in range(3) for c in (False, True)]
        cl.append(dict(theta=th, t=t, b1=b1, b2=b2, images=imgs))
    return R0, cl


def slit_z(R, c, s):
    """Point on base slit of class c at log-radial parameter s (s in [-a,0] physical, [0,a] reflected)."""
    return R * np.exp(s + 1j * (c['theta'] + c['b1'] * s**2 + c['b2'] * s**4))


def slit_dz(R, c, s):
    return slit_z(R, c, s) * (1 + 1j * (2 * c['b1'] * s + 4 * c['b2'] * s**3))


def apply_image(z, img):
    rot, conj = img
    if conj:
        z = np.conj(z)
    return z * np.exp(1j * rot)


def all_slits(R, classes):
    """list of (class index, image) for every physical slit"""
    return [(ic, img) for ic, c in enumerate(classes) for img in c['images']]
