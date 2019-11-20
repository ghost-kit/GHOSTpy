import numpy as np


def rad_to_deg(rad):
    return rad * 180/np.pi


def deg_to_rad(deg):
    return deg * np.pi/180


# Based on NASA volumetric mean radius
# http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
def cm_to_re(cm):
    return cm / 6.371008e8

# Based on NASA volumetric mean radius
# http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
def m_to_re(m):
    return m / 6.371008e6

# Based on NASA volumetric mean radius
# http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
def km_to_re(km):
    return km/6.371008e3

# Based on NASA volumetric mean radius
# http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
def re_to_cm(re):
    return re * 6.371008e8


def __rotate_x__(xyz, theta):
    Rx = np.array([[1, 0, 0], [0, np.cos(theta), np.sin(theta)], [0, -np.sin(theta), np.cos(theta)]])
    return np.dot(Rx, xyz)


def __rotate_y__(xyz, theta):
    Ry = np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
    return np.dot(Ry, xyz)


def __rotate_z__(xyz, theta):
    Rz = np.array([[np.cos(theta), np.sin(theta), 0], [-np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    return np.dot(Rz, xyz)


def __xyz_to_fan_angle__(xyz):
    pos_rot = __rotate_y__(xyz, theta=-np.pi / 2)
    x = np.arctan2(pos_rot[1], pos_rot[0])
    if x < 0:
        x += 2 * np.pi

    return x


def cart_to_polar(xy):
    r = np.sqrt(xy[0] ** 2 + xy[1] ** 2)
    theta = np.arctan2(xy[1], xy[0])
    if theta < 0:
        theta += 2 * np.pi

    return np.array([r, theta])
