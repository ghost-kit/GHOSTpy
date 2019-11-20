import numpy as np

try:
    from ghostpy.tsyganenko.t96 import t96_01
except ImportError:
    print ("T96 Not Found")
    t96_01 = None


def mon_inc(L):
    return np.all(np.diff(L) >= 0)


def mon_dec(L):
    return np.all(np.diff(L) <= 0)


def mag_s(a, axis=0):
    # print("a={} with shape: {}".format(a, np.shape(a)))
    if a is None:
        return None
    if(np.shape(a)[0] != 3):
        print ("a: {}".format(a))
        if np.shape(a)[0] != 3:
            return
    return np.sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2)


def mag(v):
    assert v is not None
    mag_shape = np.shape(v)
    if len(mag_shape) == 1:
        return mag_s(v)
    else:
        # print ("Requested: {}".format(v))
        if v.size == 0:
            return np.array([np.NaN])
        return np.apply_along_axis(mag_s, axis=1, arr=v)


def localtime_to_location(local_time):
    """
    calculates a staring location based on circle and local time
    :param local_time: Local Time to find unit vector
    :return: unit vector for seeded local time
    """
    clockangle = clock_angle(local_time)
    # need to reverse the x coordinate because... Earth.
    return np.array([-np.cos(clockangle), -np.sin(clockangle), 0.0])


def cart_to_sphere(location):
    if location is None:
        return None
    # print("location: {}".format(location))
    call_shape = np.shape(location)
    if len(call_shape) == 1:
        x = location[0]
        y = location[1]
        z = location[2]
        return cart_to_sphere_comp(x, y, z)
    else:
        return np.apply_along_axis(cart_to_sphere_p, axis=1, arr=location)


def cart_to_sphere_p(location):
    r = mag(location)
    x = location[0]
    y = location[1]
    z = location[2]
    theta = (np.pi/2) - np.arccos(z/r)
    phi = np.arctan2(y,x)
    if phi < 0 and not np.isclose(phi, 0, atol=1e-12, rtol=0):
        phi += 2*np.pi

    return np.array([r, theta, phi])


def cart_to_sphere_comp(x, y, z):
    r = mag([x,y,z])
    theta = (np.pi/2) - np.arccos(z/r)
    phi = np.arctan2(y,x)
    if phi < 0 and not np.isclose(phi, 0, atol=1e-12, rtol=0):
        phi += 2*np.pi

    return r, theta, phi


def sphere_to_cart(r, lam, phi):
    # print (np)
    x = r * np.cos(lam) * np.cos(phi)
    y = r * np.cos(lam) * np.sin(phi)
    z = r * np.sin(lam)
    return x, y, z


def clock_angle(localtime):
    """
    Generates a clock angle (in radians) from a provided local time (float)
    :param localtime: (float) representing hours on a 24 hour clock (0 - 24 ar valid)
    :return: clock angle in radians
    """
    degs = 0.25 * 60 * localtime
    return degs * np.pi / 180



def lin_interp(v, X, F, adjust=1.0):
    """
    Returns the value associated with F(v) when v is bounded by the range of X
    :param v: function parameter while looking for F(v) (v is bounded by the range of X)
    :param X: values associated with function F(X) -- values associated index by index with F
    :param F: values that define the model for F(X) -- values associated index by index with X
    :return:  linearly interpolated value for F(v) if v is in the range of X, None otherwise
    """
    # check to see if the value requested is on a node
    eq = np.where(v == X)
    if np.shape(eq)[1] >= 1:
        return F[eq][0]

    # if not, we need to interpolate
    else:
        # find the bounding points
        gt = np.where(v < X)
        lt = np.where(v > X)

        # Check to see if requested v is out of bounds, if so return None.
        if np.shape(gt)[1] == 0 or np.shape(lt)[1] == 0:
            return None

        # if within range, interpolate and return value
        if lt[0][-1] < gt[0][0]:
            lb_ind = lt[0][-1]
            ub_ind = gt[0][0]
        else:
            lb_ind = lt[0][0]
            ub_ind = gt[0][-1]

        w = (v - X[lb_ind]) / (X[ub_ind] - X[lb_ind]) * adjust
        return F[lb_ind] + ((F[ub_ind] - F[lb_ind]) * w)


