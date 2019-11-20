import inspect
import numpy as np

from ghostpy.prototype import inv_common as ic


class polarGrid:
    """
    polarGrid is a class for calculating an integration surface at a specified RE from Earth's Surface,
        where 1.0RE is the surface of the Earth.

    To Use the class, you must provide 3 things:
        1) RE: The altitude above Earth, in RE, for where the integration surface grid is to be
            initialized
        2) resolution: an array-like object in the form [longitude, latitude]
            this signifies the algorithm to generate a grid with longitude number of divisions in
            the logitude direction, and a latitude number of division in the latitude direction, starting
            at the equator (0.0) and ending at the north pole (pi/2). These arguments should both be
            integers.
        3) a function of the form 'latitude = func(longitude=x)' where both latitude and longitude
            must be represented in radians.  func(longitude=x) takes x (the longitude for which we need
            the field-line footprint) and returns the field line footprint latitude FOR THE RE THAT
            WAS SPECIFIED in the RE argument.

    NOTES: the algorithm will always pass longitude as radians.  latitude must also be returned as
            radians or the algorithm WILL FAIL.

    polarGrid offers the following user methods:

        get_low_lambda()     -- returns the lowest latitude (in radians) for the current grid space

        get_spherical_grid() -- returns a numpy meshgrid by component (grid_r, grid_l, grid_p) for
                                radius, lambda (latitude), and phi(longitude) respectively.  invalid locations of
                                the grid will be marked as a NaN in the latitude grid.

        get_cartesian_grid() -- returns the grid as a numpy meshgrid by component (grid_x, grid_y, grid_z). Invalid
                                locations within the meshgrid are marked as NaN.

        get_surface_integral(value_fun, value_coord='cart') -- provides the surface integral value for the surface
                                utilizing dual triangle cell area calculation.  The two arguments to the function are
                                1) value_fun : the function that takes a 3-tuple position and returns a value
                                    for that position
                                2) value_coord = 'cart' sets the integration value function coordinates.  Options are:
                                    'cart' (default) - coordinates will be passed to value_fun as [x,y,z] (RE)
                                    'sphere' - coordinates will be passed to value_fun as [r,lambda, phi] (RE, radians)

    """
    def __init__(self, RE=1.0, resolution=None, base_lat_fun=None):

        self.RE = RE

        # check for valid arguments
        assert resolution is not None, "POLARGRID ERROR:\nGrid resolution must be specified i" \
                                       "n the form of 'resolution=[lon,lat]'\n" \
                                       "Where:\n" \
                                       "\tlon = number of divisions in longitude\n" \
                                       "\tlat = number of divisions in latitude"

        assert len(resolution) == 2,  "POLARGRID ERROR:\nThe divs argument must contain two (and only two)" \
                                      " parts in the format [phi, lambda]" \
                                      "\nWhere phi is the number of divisions in phi (longitude)" \
                                      "\nAnd lambda is the number of divisions in lambda (latitude)"

        self.divs = resolution

        assert base_lat_fun is not None, "The base latitude function cannot be 'None'.  Please pass an " \
                                         "appropriate base latitude function.\nA proper function must respond" \
                                         " to the form:\n\n" \
                                         "base_lat = fun(longitude=valid_longitude)\n\n" \
                                         "NOTES:\n\tlongitude must be specified in radians\n" \
                                         "\tlatitude must be returned in radians"

        try:
            test = base_lat_fun(longitude=1.1)
        except TypeError:
            argument = inspect.getargspec(base_lat_fun)
            assert False, "POLARGRID ERROR:\n" \
                          "Base Latitude Function (base_lat_fun(longitude)) " \
                          "does not provide the correct arguments.\n" \
                          "The function must respond to a call of function(longitude=valid_longitude)\n" \
                          "where longitude is in radians.\n\n" \
                          "The function provided has the arguments: {}\n".format(argument.args)

        self.base_lat_fun = base_lat_fun

        # define the grid structure
        self.grid_phi = np.linspace(start=0, stop=2*np.pi, num=self.divs[0], dtype=np.float64, endpoint=False)

        self.base_lambdas = [base_lat_fun(longitude=x) for x in self.grid_phi]
        assert np.all(np.array(self.base_lambdas) <= np.pi/2), "POLARGRID ERROR WHILE BUILDING POLAR CAP GRID:\n" \
                                                               "There is a problem with your starting " \
                                                               "Lambdas (Latitude). One or more are out of bounds\n" \
                                                               "Please check your base_lat_fun(phi) function to " \
                                                               "ensure that it is producing the correct results.\n" \
                                                               "Check to make sure you are returning values in radians"\
                                                               "\nResults are limited to between 0 " \
                                                               "and pi/2 radians\n\n" \
                                                               "NOTE: Floating point calculation errors can trigger " \
                                                               "this error. Correct for this in your function."

        self.grid_lambda = np.linspace(start=self.get_low_lambda(), stop=np.pi / 2, num=self.divs[1], endpoint=True, dtype=np.float64)

        self.phi_diff = self.grid_phi[1] - self.grid_phi[0]

    def get_low_lambda(self):
        """
        Method for retrieving the lowest latitude on the grid.
        :return: lowest latitude (lambda) on the grid (in radians)
        """
        return np.min(self.base_lambdas)

    def get_spherical_grid(self):
        """
        Method for retrieving the meshgrid in spherical coordinates
        :return: grid_r, grid_l, grid_p (r, lambda, phi) (i.e. radius, latitude, longitude)
        """
        base_grid_r, base_grid_l, base_grid_p = np.meshgrid(self.RE,
                                                            self.grid_lambda,
                                                            self.grid_phi,
                                                            indexing='ij')

        # Mark cells that are not within the scope as lowest lambda
        # and adjust the lowest level to be the min for the given phi
        for x in range(len(self.base_lambdas)):
            invalid_loc = np.argwhere(base_grid_l[0, :, x] < self.base_lambdas[x])
            # print invalid_loc
            if len(invalid_loc > 0):
                # base_grid_l[0, invalid_loc, x] = np.nan
                base_grid_l[0,invalid_loc, x] = self.base_lambdas[x]

        return base_grid_r, base_grid_l, base_grid_p

    def get_cartesian_grid(self):
        """
        Method for acquiring a cartesian version of the grid.
        :return: grid_x, grid_y, grid_z
        """
        r, l, p = self.get_spherical_grid()
        return ic.get_location_from_theta_phi_r(r=r, theta=l, phi=p)

    def get_surface_integral(self, value_fun=None, value_coord='cart'):
        """
        method to select between two versions of te integration methods
        :param value_fun:
        :param value_coord:
        :return:
        """
        assert value_coord.lower() == 'cart' or value_coord.lower() == 'sphere', "ERROR: Valid Value for value_coord " \
                                                                                 "are: \n\n\t'cart' \n\t'sphere"

        if value_fun is None:
            value_fun = lambda xr, yl, zp: np.ones(np.shape(xr))

        arguments = inspect.getargspec(value_fun).args

        assert np.any(np.array(arguments) == 'xr') and \
            np.any(np.array(arguments) == 'yl') and \
            np.any(np.array(arguments) == 'zp'), "ERROR: value_fun(xr, yl, zp) must accept " \
                                                 "all arguments.\nxr = [numpy array] grid x (cart) or r (sphere) compontents\n" \
                                                 "yl = [numpy array] grid y (cart) or lambda [latitude] (sphere) compontents\n" \
                                                 "zp = [numpy array] grid z (cart) or phi [longitude] (shpere) compontents"

        if value_coord.lower() == 'sphere':
            conv = True
        else:
            conv = False

        grid_r, grid_l, grid_p = self.get_spherical_grid()

        # get lambda diffs
        lambda_diffs = np.apply_along_axis(np.diff, 0, grid_l[0,:,:], n=1)

        area_bot = self.RE**2 * np.cos(grid_l[0, :-1, :]) * lambda_diffs * self.phi_diff
        area_top = self.RE**2 * np.cos(grid_l[0, 1:, :]) * lambda_diffs * self.phi_diff

        area = (0.5*area_bot) + (0.5*area_top)

        if value_coord.lower() == "cart":
            grid_1, grid_2, grid_3 = self.get_cartesian_grid()
        else:
            grid_1 = grid_r
            grid_2 = grid_l
            grid_3 = grid_p

        v1 = value_fun(xr=grid_1, yl=grid_2, zp=grid_3)

        vtop1 = v1[0,1:,:]
        vbot1 = v1[0,:-1,:]
        vtop2 = np.roll(v1,-1, axis=2)[0,1:,:]
        vbot2 = np.roll(v1,-1, axis=2)[0,:-1,:]

        val = (vtop1 + vtop2 + vbot1 + vbot2)/4

        integral = np.sum(area * val)

        return integral

