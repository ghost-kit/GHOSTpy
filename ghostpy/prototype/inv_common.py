import bisect as bs
import collections as col
import inspect
from types import *

import numpy as np
from scipy import interpolate

try:
    from ghostpy.tsyganenko.t96 import t96_01
except ImportError:
    print ("Failure getting T96")
    pass

def mon_inc(L):
    return np.all(np.diff(L) >= 0)

def mon_dec(L):
    return np.all(np.diff(L) <= 0)



def mag(vector):
    tot = np.array(vector, dtype=np.float64)**2
    return np.sqrt(np.sum(tot, dtype=np.float64))


def dict_to_x_y(dict_to_plot={}):
    x = sorted(dict_to_plot.keys())
    y = [dict_to_plot[k] for k in x]
    return x, y


def data_to_dict(**kwargs):
    x = None
    y = None
    z = None
    new_dict = col.OrderedDict()

    for key, value in kwargs.iteritems():
        nv = np.array(value, dtype=np.float64)
        try:
            if len(nv[0]) == 3:
                x = [d[0] for d in nv]
                y = [d[1] for d in nv]
                z = [d[2] for d in nv]
                key_x = str(key) + str("_X")
                key_y = str(key) + str("_Y")
                key_z = str(key) + str("_Z")
                new_dict[key_x] = x
                new_dict[key_y] = y
                new_dict[key_z] = z
            else:
                new_dict[str(key)] = np.array(nv, dtype=np.float64)
        except:
            new_dict[str(key)] = np.array(nv, dtype=np.float64)

    return new_dict


def lin_interp(A, B, v):
    """
    Do a linear Interpolation between points on A. v is a value that exists between the first and last element of B.
    :param A: Points that are defined. We want to get a point in between two of these values.
    :param B: Values associated with Points in B.
    :param v: Value that we are looking for the point of.
    :return: Point found from the interpolation
    """
    # print "running linear interpolation"
    if mon_inc(B) or mon_dec(B):
        # print "Array B is good to go!"
        if in_range(B, v):
            # print "Value is in range!"
            if np.isclose(B[0], v, atol=1e-8, rtol=0.0):
                # print "close enough"
                return A[0]
            else:
                hi_ind = bs.bisect_left(B, v)
                lo_ind = hi_ind -1
                W = (v - B[lo_ind])/(B[hi_ind] - B[lo_ind])

                # print "hi index: {}\nlow index: {}\nw: {}".format(hi_ind, lo_ind, W)

                return np.array(A[lo_ind]) + (np.array((A[hi_ind]) - np.array(A[lo_ind])) * W)
        else:
            print ("Value is out of range")
            print ("Value given: {}".format(v))
            print ("Range: {} - {}".format(B[0], B[-1]))
            return None
    else:
        print ("What happened?")
        print ("Mon-inc? {}".format(mon_inc(B)))
        print ("Mon-dec? {}".format(mon_dec(B)))
        print (B)
        return None


def cub_interp(A, B, v):
    """
    Qucik convenience method to wrap cubic spline to work same as the linear interpolator above
    :param A:
    :param B:
    :param v:
    :return:
    """
    f = interpolate.splrep(B, A, s=0)
    return interpolate.splev(v, f, der=0)


def in_range(A, v):
    """
    Checks to see if v is in the RANGE of A.
    :param A: Vector-Like object (of scalars)
    :param v: scalar value
    :return: True if v is in Range of A, else False
    """
    return A[0] <= v <= A[-1]


def get_location_from_localtime(local_time):
    """
    calculates a staring location based on circle and local time
    :param local_time: Local Time to find unit vector
    :return: unit vector for seeded local time
    """
    clockangle = _get_clock_angle(local_time)
    # need to reverse the x coordinate because... Earth.
    return np.array([-np.cos(clockangle), -np.sin(clockangle), 0.0])


def get_r_theta_phi_from_location(location):
    x = location[0]
    y = location[1]
    z = location[2]
    return get_r_theta_phi_from_x_y_z(x,y,z)


def get_r_theta_phi_from_x_y_z(x,y,z):
    r = mag([x,y,z])
    theta = (np.pi/2) - np.arccos(z/r)
    phi = np.arctan2(y,x)

    if phi < 0:
        phi += 2*np.pi

    return r, theta, phi


def get_location_from_theta_phi_r(theta, phi, r):
    x = r * np.cos(theta) * np.cos(phi)
    y = r * np.cos(theta) * np.sin(phi)
    z = r * np.sin(theta)

    # if np.isnan(x) or np.isnan(y) or np.isnan(z):
    #     print "x: {}, y: {}, z:{}".format(x,y,z)
    #     print "r: {}, theta: {}, phi: {}".format(r, theta, phi)
    #     sys.exit()

    return x, y, z




def _get_clock_angle(localtime):
    """
    Generates a clock angle (in radians) from a provided local time (float)
    :param localtime: (float) representing hours on a 24 hour clock (0 - 24 ar valid)
    :return: clock angle in radians
    """
    degs = 0.25 * 60 * localtime
    return degs * np.pi / 180


def get_local_time_from_location(location):
    """
    returns the local time from a given location (location != [0,0,z])
    :param location: [x,y,z] coordinates to translate. Cannot be [0,0,z]
    :return: local time as a float
    """
    RE = mag(location)
    unit = np.array(location)/RE
    clockangle = np.arctan(unit[1]/unit[0])
    ca_deg = (clockangle * 180/np.pi)

    # Because of how ArcTan works, we need to adjust the calculated angle based on what quadrant it is in.
    if location[0] < 0 and location[1] <= 0:
        pass
    elif location[0] < 0 and location[1] >= 0:
        ca_deg += +360
    elif location[0] >=0 and location[1] < 0:
        ca_deg += 180
    else:
        ca_deg += 180

    time = ca_deg/(0.25*60)

    return time


def lon_to_lt(lon):
    return phi_to_lt(deg_to_rad(lon))


def phi_to_lt(phi):
    loc = get_location_from_theta_phi_r(r=1.0, theta=0, phi=phi)
    return get_local_time_from_location(loc)


def lt_to_lon(lt):
    loc = get_location_from_localtime(lt)
    return get_lon(loc)


def get_lat(xyz):
    rtp = get_r_theta_phi_from_location(xyz)
    return rad_to_deg(rtp[1])


def get_lon(xyz):
    rtp = get_r_theta_phi_from_location(xyz)
    return rad_to_deg(rtp[2])


def rad_to_deg(rad):
    return rad * 180/np.pi


def deg_to_rad(deg):
    return deg * np.pi/180


def build_spherical_cap_grid(RE=1.0, base_lat=0, divs=100, grid_type="sphere"):
    # TODO: REMOVE THIS FUNCTION ONCE 2 (below) is online
    """
    Generate a spherical cap grid with cells = divs^2
    :param RE:
    :param base_lat:
    :param divs:
    :param grid_type:
    :return:
    """
    lons = np.linspace(start=0, stop=360, num=divs, endpoint=False)
    lats = np.linspace(start=base_lat, stop=90, num=divs, endpoint=True)

    cart_grid = col.OrderedDict()
    sphere_grid = col.OrderedDict()

    for lon in lons:
        lt = lon_to_lt(lon)
        cart_grid[lt] = []
        sphere_grid[lt] = []

        for base_lat in lats:
            theta = deg_to_rad(base_lat)
            phi = deg_to_rad(lon)
            if np.isnan(phi) or np.isnan(theta):
                print "theta: {}, phi: {}".format(theta, phi)
                sys.exit()
            # print theta, phi
            # print "pi/2: {}".format(np.pi/2)
            if theta > np.pi/2:
                theta = np.pi/2
                sys.exit()
            pointc = get_location_from_theta_phi_r(theta=theta, phi=phi, r=RE)
            points = [RE, theta, phi]
            # print points
            if grid_type is 'cart':
                cart_grid[lt].append(pointc)
            elif grid_type is 'sphere':
                sphere_grid[lt].append(points)
            elif grid_type is 'both':
                cart_grid[lt].append(pointc)
                sphere_grid[lt].append(points)
            else:
                print "{} not an understood grid type".format(grid_type)
                sys.exit()

    if grid_type is 'cart':
        return cart_grid
    elif grid_type is 'sphere':
        return sphere_grid
    elif grid_type is 'both':
        return cart_grid, sphere_grid


def build_spherical_cap_grid2(RE=1.0, base_lat_fun=None, divs=None, grid_type="sphere"):
    # This is a refactoring of the above "build_spherical_cap_grid" method to utilize
    #   a base latitude function.
    """
    Method for the construction of a polar cap grid at a given RE distance from the center of the Earth. RE=1.0
    represents the surface of the earth.

    Example usage:
        X, Y, Z = build_spherical_cap_grid2(RE=2.5, base_lat_fun=lat_lookup, divs = [100,50], grid_type='cart')

    This will return a Numpy meshgrid for X, Y, and Z that make up a polar cap grid that starts at
        base_lat_fun(phi) for each phi in the grid, and will have a base resolution of divs[0] x divs[1]
        where divs[0] are the number of phi divisions, and divs[1] represents the  number of lambda divisions
        for the lowest lambda in the set. All lambda divsions are at the same latitude for each phi, with the
        exception of the base latitude, which is defined by the base_lat_fun(phi) function.

    :param RE: Number of Earth Radii that defines the radius of the spherical cap

    :param base_lat_fun: Function in the form lat_fun(localtime) to determine the latitude for the
            field line intersection point at the given local time.

    :param divs: the base number of latitude and longitude divisions for the grid.  The actual
            number of divisions is determined by the lambda of the given localtime. The full
            number of divisions will be utilized on the lowest lambda intersection.

    :param grid_type: type of grid to return. "sphere" for spherical, "cart" for cartesian.

    :return:  Returns the type of grid specified, with the lower boundaries defined by the
            base_lat_fun function, and the resolution defined by "divs"
    """

    assert type(grid_type) is StringType, "ERROR:\nGrid type argument must be a string.\n" \
                                          "The options are:\n" \
                                          "\t'sphere' for a grid in spherical coordinates\n" \
                                          "\t'cart' for a grid in Cartesian coordinates"

    if divs is None:
        divs = np.array([100, 100])
    else:
        # check for valid divisions
        assert len(divs) == 2, "ERROR:\nThe divs argument must contain two (and only two)" \
                               " parts in the format [phi, lambda]" \
                               "\nWhere phi is the number of divisions in phi (longitude)" \
                               "\nAnd lambda is the number of divisions in lambda (latitude)"

    # check to make sure a function is defined
    assert base_lat_fun is not None
    try:
        test = base_lat_fun(phi=0.0)
    except:
        print inspect.getargspec(base_lat_fun)
        assert False, "ERROR: Base Latitude Function (base_lat_fun(phi))" \
                      " does not provide the correct arguments.\n" \
                      "The function must respond to a call of function(phi=longitude)\n" \
                      "where phi is in radians"

    grid_phi = np.linspace(start=0, stop=2*np.pi, num=divs[0])

    start_lambdas = [base_lat_fun(phi=x) for x in grid_phi]
    assert np.all(np.array(start_lambdas) <= np.pi/2), "ERROR WHILE BUILDING POLAR CAP GRID:\n" \
                                                       "There is a problem with your starting Lambdas (Latitude). " \
                                                       "One or more are out of bounds\n" \
                                                       "Please check your base_lat_fun(phi) function to " \
                                                       "ensure that it is producing the correct results.\n" \
                                                       "Results are limited to between 0 and pi/2"

    low_lambda = np.min(start_lambdas)
    grid_lambda = np.linspace(start=low_lambda, stop=np.pi/2, num=divs[1], endpoint=True)
    grid_r, grid_l, grid_p = np.meshgrid(RE, grid_lambda, grid_phi, indexing='ij')

    # Mark cells that are not within the scope as NaN
    # and adjust the lowest level to be the min for the given phi
    for x in range(len(start_lambdas)):
        invalid_loc = np.argwhere(grid_l[0,:,x] < start_lambdas[x])
        # print invalid_loc
        if len(invalid_loc > 0):
            min_loc = invalid_loc[-1]
            grid_l[0,min_loc,x] = start_lambdas[x]
            grid_l[0, invalid_loc[:-1], x] = np.nan

    if grid_type.lower() == "sphere":
        return grid_r, grid_l, grid_p
    elif grid_type.lower() == "cart":
        return get_location_from_theta_phi_r(r=grid_r, theta=grid_l, phi=grid_p)
    else:
        print "type unknown"
        return None, None, None


def integrate_flux_cap_area(sphere_grid=None, pvalue_fun=None):
    # TODO: integrate this integration into driftShell calculations. It is a cleaner function.
    """
    integrate flux over the given grid (rectangular integration)
    :param grid: spherical cap grid
    :param data_fun: function that returns data at the r,phi,theta on the grid
    :return: value of the integration
    """

    if pvalue_fun is None:
        pvalue_fun = lambda x1, x2, x3, x4: 1

    if sphere_grid is None:
        print "a grid is required to integrate"
        sys.exit()

    summation = 0
    elements = {}

    hkeys = sorted(sphere_grid.keys())
    lons = [lt_to_lon(a) for a in hkeys]
    phis = [deg_to_rad(a) for a in lons]
    for key in hkeys:
        elements[key] = []

    num_lon = len(hkeys)

    for index in range(len(hkeys)):
        # need to integrate by phi
        if index+1 < len(hkeys):
            i1 = index
            i2 = i1 + 1
        else:
            i1 = index
            i2 = 0

        lt_lats1 = sphere_grid[hkeys[i1]]
        lt_lats2 = sphere_grid[hkeys[i2]]

        for index2 in range(len(lt_lats1)-1):
            point_lon1_lat1 = lt_lats1[index2]
            point_lon1_lat2 = lt_lats1[index2+1]
            point_lon2_lat1 = lt_lats2[index2]
            point_lon2_lat2 = lt_lats2[index2+1]

            t_ln1_lt1 = point_lon1_lat1[1]
            t_ln1_lt2 = point_lon1_lat2[1]
            t_ln2_lt1 = point_lon2_lat1[1]
            t_ln2_lt2 = point_lon2_lat2[1]

            p_ln1_lt1 = point_lon1_lat1[2]
            p_ln1_lt2 = point_lon1_lat2[2]
            p_ln2_lt1 = point_lon2_lat1[2]
            p_ln2_lt2 = point_lon2_lat2[2]

            r = point_lon1_lat1[0]

            dT1 = t_ln1_lt2 - t_ln1_lt1

            if p_ln2_lt1 < p_ln1_lt1:
                # print p_ln2_lt1 - p_ln1_lt1
                dP = (p_ln2_lt1 + (2 * np.pi)) - p_ln1_lt1
                # print dP
            else:
                dP = p_ln2_lt1 - p_ln1_lt1

            # print "dP: {}".format(dP)
            # print "dT: {}".format(dT)
            # print "R: {}".format(r)
            area = r**2 * np.cos(point_lon1_lat1[1]) * dT1 * dP
            if area < 0:
                print "invalid area: {}".format(area)
                sys.exit()
            # print "area: {}".format(area)
            elements[hkeys[index]].append(area * pvalue_fun(point_lon1_lat1, point_lon1_lat2, point_lon2_lat1, point_lon2_lat2))

    for key in elements:
        summation += np.sum(elements[key])

    return summation

def integrate_flux_cap_area2(sphere_grid=None, pvalue_fun=None, error_tol=1e-8):
    """
    integrate flux over the given grid (Triangular integration)
    :param grid: spherical cap grid
    :param data_fun: function that returns data at the r,phi,theta on the grid
    :return: value of the integration
    """

    if pvalue_fun is None:
        pvalue_fun = lambda x1, x2, x3, x4: 1

    if sphere_grid is None:
        print "a grid is required to integrate"
        sys.exit()

    summation = 0
    elements = {}

    hkeys = sorted(sphere_grid.keys())
    lons = [lt_to_lon(a) for a in hkeys]
    phis = [deg_to_rad(a) for a in lons]
    for key in hkeys:
        elements[key] = []

    num_lon = len(hkeys)

    for index in range(len(hkeys)):
        # need to integrate by phi
        if index+1 < len(hkeys):
            i1 = index
            i2 = i1 + 1
        else:
            i1 = index
            i2 = 0

        lt_lats1 = sphere_grid[hkeys[i1]]
        lt_lats2 = sphere_grid[hkeys[i2]]

        for index2 in range(len(lt_lats1)-1):
            point_lon1_lat1 = lt_lats1[index2]
            point_lon1_lat2 = lt_lats1[index2+1]
            point_lon2_lat1 = lt_lats2[index2]
            point_lon2_lat2 = lt_lats2[index2+1]

            t_ln1_lt1 = point_lon1_lat1[1]
            t_ln1_lt2 = point_lon1_lat2[1]
            t_ln2_lt1 = point_lon2_lat1[1]
            t_ln2_lt2 = point_lon2_lat2[1]

            p_ln1_lt1 = point_lon1_lat1[2]
            p_ln1_lt2 = point_lon1_lat2[2]
            p_ln2_lt1 = point_lon2_lat1[2]
            p_ln2_lt2 = point_lon2_lat2[2]

            r = point_lon1_lat1[0]

            dT1 = t_ln1_lt2 - t_ln1_lt1
            dT2 = t_ln2_lt2 - t_ln2_lt1

            if p_ln2_lt1 < p_ln1_lt1:
                # print p_ln2_lt1 - p_ln1_lt1
                dP1 = (p_ln2_lt1 + (2 * np.pi)) - p_ln1_lt1
                # print dP
            else:
                dP1 = p_ln2_lt1 - p_ln1_lt1

            if p_ln2_lt2 < p_ln1_lt2:
                dP2 = (p_ln2_lt2 + (2 * np.pi)) - p_ln1_lt2
            else:
                dP2 = p_ln2_lt2 - p_ln1_lt2

            area1 = 0.5*(r**2 * np.cos(point_lon1_lat1[1]) * dT1 * dP1)
            if area1 < 0:
                if area1 < 0:
                    if np.isclose(0, area1, rtol=0, atol=error_tol):
                        area1 = 0.0
                    else:
                        print "invalid area: {}".format(area1)
                        sys.exit()
            area2 = 0.5*(r**2 * np.cos(point_lon1_lat2[1]) * dT2 * dP2)
            if area2 < 0:
                if np.isclose(0, area2, rtol=0, atol=error_tol):
                    area2 = 0.0
                else:
                    print "invalid area: {}".format(area2)
                    sys.exit()
            elements[hkeys[index]].append((area1 + area2) * pvalue_fun(point_lon1_lat1, point_lon1_lat2, point_lon2_lat1, point_lon2_lat2))
    for key in elements:
        summation += np.sum(elements[key])

    return summation


def analytic_sphere_cap_area(lat, r):
    """
    Calculate the area of the spherical cap based on latitude and R
    :param lat:
    :param r:
    :return:
    """
    H = r - np.sin(deg_to_rad(lat))
    return 2 * np.pi * r * H


def analytic_dipole_L(lat, r):
    theta = deg_to_rad(lat)
    return r/np.cos(theta)**2


def analytic_dipole_line_lat_crossing(L, r):
    Ls = np.array(L)
    sin_theta = np.sqrt(r/Ls)
    theta = np.arccos(sin_theta)

    return rad_to_deg(theta)


def get_dipole_value(XYZ, dipole_tilt=0.0, Be=3.15e4):
    X = np.array([A[0] for A in XYZ])
    Y = np.array([A[1] for A in XYZ])
    Z = np.array([A[2] for A in XYZ])

    P = X**2
    U = Z**2
    T = Y**2
    Q = Be/np.sqrt(P + U + T)**5
    V = 3.0 * Z * X
    sps = np.sin(dipole_tilt)
    cps = np.cos(dipole_tilt)

    Bx = Q * ((T + U - 2.0 * P) * sps - V * cps)
    By = -3.0*Y * Q * (X * sps + Z * cps)
    Bz = Q * ((P + T - 2.0 * U) * cps - V * sps)

    tup = []
    for ind in range(len(Bx)):
        tup.append((Bx[ind], By[ind], Bz[ind]))

    return tup


def trace_rk45(inner_boundary=0.5, h=1e-3, x0=None, max_steps=2000, error_tol=1e-4, direction="f"):
    if x0 is None:
        x0 = [6.0, 0.0, 0.0]

    mag_x0 = mag(x0)
    x = x0
    RE = mag_x0
    dpv = dpval(x)
    path = [tuple(x)]
    value = [tuple(dpv)]
    steps = 0
    while RE > inner_boundary:
        s = 0
        while not np.isclose(s, 1.0, atol=0.0001, rtol=0.00) and not np.isnan(s):
            k1 = h * dpval(x, direction=direction)
            k2 = h * dpval(x + 1. / 4. * k1, direction=direction)
            k3 = h * dpval(x + 3. / 32. * k1 + 9. / 32. * k2, direction=direction)
            k4 = h * dpval(x + 1932. / 2197. * k1 - 7200. / 2197. * k2 + 7296. / 2197. * k3, direction=direction)
            k5 = h * dpval(x + 439. / 216. * k1 - 8. * k2 + 3680. / 513. * k3 - 845. / 4104. * k4, direction=direction)
            k6 = h * dpval(x - 8. / 27. * k1 + 2. * k2 - 3544. / 2565. * k3 + 1859. / 4104. * k4 - 11. / 40. * k5, direction=direction)

            x1 = x + 25. / 216. * k1 + 1408. / 2565. * k3 + 2197. / 4104. * k4 - 1. / 5. * k5
            x2 = x + 16. / 135. * k1 + 6656. / 12825. * k3 + 28561. / 56430. * k4 - 9. / 50. * k5 + 2. / 55. * k6

            # print "x2: {} - x1: {} = {}".format(x2,x1,x2-x1)

            s = 0.84 * (error_tol * h / (mag(x2 - x1))) ** 0.25

            if np.isnan(s) or np.isinf(s):
                s = 1.0
            h *= s
            # print "s: {}".format(s)

        x = x1
        RE = mag(x)
        dpv = dpval(x)

        # print "X: {}  ({}RE)... h = {}".format(x, RE, h)
        steps += 1

        if steps > max_steps:
            print "Truncating on number of steps.\nConsider Increasing number of steps."
            break
        path.append(tuple(x1))
        value.append(tuple(dpv))
        # print "location: {} (at {} RE) with Step Size: {}".format(x, dpv, RE, h)

    return path, value


def dpval(x, direction='f'):
    data = np.array(get_dipole_value([x])[0])
    if direction is 'f':
        return data
    elif direction is 'b':
        return -data


def dpval2(x, direction='f'):
    data = np.array(get_dipole_value2(x))
    if direction is 'f':
        return data
    elif direction is 'b':
        return -data


def cell_avg_dipole(x1, x2, x3, x4):
    point1 = np.array(get_location_from_theta_phi_r(r=x1[0], theta=x1[1], phi=x1[2]))
    point2 = np.array(get_location_from_theta_phi_r(r=x2[0], theta=x2[1], phi=x2[2]))
    point3 = np.array(get_location_from_theta_phi_r(r=x3[0], theta=x3[1], phi=x3[2]))
    point4 = np.array(get_location_from_theta_phi_r(r=x4[0], theta=x4[1], phi=x4[2]))

    pointC = (point1 + point2 + point3 + point4)/4
    pointC /= mag(pointC)

    data1 = np.array(dpval(point1))
    data2 = np.array(dpval(point2))
    data3 = np.array(dpval(point3))
    data4 = np.array(dpval(point4))

    dataC = np.array(dpval(pointC))

    dp1 = np.dot(data1, point1)
    dp2 = np.dot(data2, point2)
    dp3 = np.dot(data3, point3)
    dp4 = np.dot(data4, point4)
    dpC = np.dot(dataC, pointC)

    dp_avg = (dp1 + dp2 + dp3 + dp4 + dpC)/5
    # dp_avg = (dp1 + dp2 + dp3 + dp4)/4
    return dp_avg


def cell_cntr_dipole(x1, x2, x3, x4):
    point1 = np.array(get_location_from_theta_phi_r(r=x1[0], theta=x1[1], phi=x1[2]))
    point2 = np.array(get_location_from_theta_phi_r(r=x2[0], theta=x2[1], phi=x2[2]))
    point3 = np.array(get_location_from_theta_phi_r(r=x3[0], theta=x3[1], phi=x3[2]))
    point4 = np.array(get_location_from_theta_phi_r(r=x4[0], theta=x4[1], phi=x4[2]))

    pointC = (point1 + point2 + point3 + point4)/4
    pointC /= mag(pointC)

    # data1 = np.array(dpval(point1))
    # data2 = np.array(dpval(point2))
    # data3 = np.array(dpval(point3))
    # data4 = np.array(dpval(point4))

    dataC = np.array(dpval(pointC))

    # dp1 = np.dot(data1, point1)
    # dp2 = np.dot(data2, point2)
    # dp3 = np.dot(data3, point3)
    # dp4 = np.dot(data4, point4)
    dpC = np.dot(dataC, pointC)

    dp_avg = dpC

    return dp_avg


def triangle_dipole(x1, x2, x3):
    point1 = np.array(get_location_from_theta_phi_r(r=x1[0], theta=x1[1], phi=x1[2]))
    point2 = np.array(get_location_from_theta_phi_r(r=x2[0], theta=x2[1], phi=x2[2]))
    point3 = np.array(get_location_from_theta_phi_r(r=x3[0], theta=x3[1], phi=x3[2]))

    pointC = (point1 + point2 + point3)/3
    pointC /= mag(pointC)

    point12 = (point1 + point2)/2
    point23 = (point2 + point3)/2
    point13 = (point1 + point3)/2

    point12 /= mag(point12)
    point23 /= mag(point23)
    point13 /= mag(point13)


    data1 = np.array(dpval(point1))
    data2 = np.array(dpval(point2))
    data3 = np.array(dpval(point3))
    dataC = np.array(dpval(pointC))

    data12 = np.array(dpval(point12))
    data13 = np.array(dpval(point13))
    data23 = np.array(dpval(point23))

    dp1 = np.dot(data1, point1)
    dp2 = np.dot(data2, point2)
    dp3 = np.dot(data3, point3)
    dpC = np.dot(dataC, pointC)
    d12 = np.dot(data12, point12)
    d13 = np.dot(data13, point13)
    d23 = np.dot(data23, point23)

    dp_avg = (dp1 + d12 + (2*d13) + d23 + (3*dpC) + dp2 + dp3)/10

    return dp_avg


def L_star_from_flux(flux, Be=-3.15e4, RE=1.0):
    Rearth = RE
    L = 2 * np.pi * Be * Rearth / flux
    return L


def get_dipole_value2(XYZ, r=1, Be=-3.15e4):
    r5 = r**-5
    bx = 3 * XYZ[0] * XYZ[2] * Be * r5
    by = 3 * XYZ[1] * XYZ[2] * Be * r5
    bz = (3 * XYZ[2]**2 - r**2) * Be * r5

    return bx, by, bz


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


def t96_field(x, y, z, p_dyn, dst, imf_by, imf_bz, ps):
    """
    Generates the T96 model fields and modifies a dipole
    :param p_dyn: Solar Wind Dynamic Pressure
    :param dst: Disturbance (storm time) index.
    :param imf_by: Interplanetary Magnetic Field - Y component (GSM)
    :param imf_bz: Interplanetary Magnetic Field - Z component (GSM)
    :param ps:
    :return:
    """
    # get the dipole
    b_x, b_y, b_z = dipole_field(x=x, y=y, z=z, ps=ps)

    # add the t96 external components
    print "Dst: ", dst
    print "Pdyn: ", p_dyn
    print "IMF_By: ", imf_by
    print "IMF_Bz: ", imf_bz
    print "dipole_tilt: ", ps

    bx, by, bz = t96_wrap(x=x, y=y, z=z, p_dyn=p_dyn, dst=dst, imf_by=imf_by, imf_bz=imf_bz, ps=ps)
    b_x = b_x + np.array(bx)
    b_y = b_y + np.array(by)
    b_z = b_z + np.array(bz)

    return b_x, b_y, b_z


def stokes_integration_r_theta_phi(path, val_fun=None, Be=-3.15e4):
    if val_fun is None:
        val_fun = lambda x, y, z: 1

    summation = 0

    for index in range(len(path)):
        if index + 1 >= len(path):
            ind0 = index
            ind1 = 0
        else:
            ind0 = index
            ind1 = index + 1

        ds2 = (path[ind1][2] - path[ind0][2]) * path[ind0][0]
        # fix for loop back
        if ds2 < 0:
            ds2 += (np.pi *2)

        summation += Be * path[ind0][0]**2 * (np.sin((path[ind1][1]+path[ind0][1])/2)**2) * ds2

    return summation


def t96_wrap(x, y, z, p_dyn, dst, imf_by, imf_bz, ps):
    parmod = np.array([p_dyn, dst, imf_by, imf_bz, 0, 0, 0, 0, 0, 0])
    print "Calculating the T96 Model... please stand by..."
    t96_v = np.vectorize(t96_01, excluded=['parmod', 'ps'])

    return t96_v(x=x, y=y, z=z, parmod=parmod, ps=ps)


def get_lat_path(lat=0, divs=24, re=1):
    phi_list = np.linspace(start=0, stop=np.pi*2, num=divs, endpoint=False)
    path_list = []
    theta = (np.pi/2)-deg_to_rad(lat)
    for phi in phi_list:
        path_list.append((re, theta, phi))

    return path_list

def re_to_cm(re):
    return re * 6.38e8

def cm_to_re(cm):
    return cm / 6.38e8

def m_to_re(m):
    return m / 6.38e6

def km_to_re(km):
    return km / 6.38e3
