import numpy as np
import ghostpy.algorithms.common as algc


def newell_surf_normal(poly):
    # This algorithm is  from the Newell's method pseudo-code found at:
    # https://www.khronos.org/opengl/wiki/Calculating_a_Surface_Normal

    old_settings = np.geterr()
    np.seterr(invalid='ignore')
    norm = np.array([0.0, 0.0, 0.0])
    for i in range(len(poly)):
        v_curr = poly[i]
        v_next = poly[(i + 1) % len(poly)]
        norm[0] += (v_curr[1] - v_next[1]) * (v_curr[2] + v_next[2])
        norm[1] += (v_curr[2] - v_next[2]) * (v_curr[0] + v_next[0])
        norm[2] += (v_curr[0] - v_next[0]) * (v_curr[1] + v_next[1])

    mag_norm = algc.mag(norm)
    unit_norm = norm / mag_norm

    np.seterr(invalid=old_settings['invalid'])
    return unit_norm


def cross_surf_norm(p1, p2, p3):
    l1 = vector_between(p1, p2)
    l2 = vector_between(p1, p3)
    cp = np.cross(l1, l2)
    mcp = algc.mag(cp)
    return cp / mcp


def perp_dist(unit_normal, plane_point, non_plane_point):
    """
    Provides a vector from a point to a plane that is perpendicular to the plane
    :param unit_normal: unit normal vector of the plane
    :param plane_point: any point on the plane
    :param non_plane_point: the point you want the vector from
    :return: the perpendicular vector from non_plane_point to the surface.
    """
    f = vector_between(plane_point, non_plane_point)
    dist = np.dot(f, unit_normal)
    print ("dist: {}".format(dist))
    return -dist


def vector_between(a, b):
    return b-a


def vector_direction(a, b):
    return np.dot(a,b)
