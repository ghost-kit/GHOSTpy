import numpy as np
import ghostpy.algorithms.common as algc
import ghostpy.algorithms.convert as algx


def dipole_field(x, y, z, ps=0., Be=3.15e4):
    """
    Generates a dipole field on the internal grid.
    :param ps: Dipole Tilt. Default is 0.
    :return: X, Y, Z components of the dipole field for the internal grid (GSM Coordinates)
    """
    P = x ** 2
    U = z ** 2
    T = y ** 2
    Q = Be / np.sqrt(P + U + T) ** 5
    V = 3.0 * z * x
    sps = np.sin(ps)
    cps = np.cos(ps)

    Bx = Q * ((T + U - 2.0 * P) * sps - V * cps)
    By = -3.0 * y * Q * (x * sps + z * cps)
    Bz = Q * ((P + T - 2.0 * U) * cps - V * sps)
    return Bx, By, Bz

def dipole_L(lat, r):
    theta = algx.deg_to_rad(lat)
    return r/np.cos(theta)**2

def dipole_L_rad(lam,r):
    return r/np.cos(lam)**2



def dipole_footprint(L, r):
    Ls = np.array(L)
    sin_theta = np.sqrt(r/Ls)
    theta = np.arccos(sin_theta)

    return theta