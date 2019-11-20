from __future__ import absolute_import

import numpy as np
import unittest as ut

from grids import polarGrid as pg
from prototype import inv_common as ic

# These are timing variables
setup_string = 'from grids import polarGrid as pg\n' \
               'import numpy as np'
# num_samples = 100
res = [75,150]
RE=2.5

def test_fun(longitude):
    assert longitude is not None
    # return longitude/4.0
    return 0

def test_value_fun(xr, yl, zp):
    dpX, dpY, dpZ = ic.dipole_field(xr, yl, zp)

    b = np.vstack(([dpX.T], [dpY.T], [dpZ.T])).T
    n = np.vstack(([xr.T], [yl.T], [zp.T])).T

    # print("B Values: {}".format(b))
    # print("N Values: {}".format(n))

    val = np.zeros_like(dpX)
    shape = np.shape(dpX)
    for x in range(shape[0]):
        for y in range(shape[1]):
            for z in range(shape[2]):
                val[x,y,z] = b[x,y,z].dot(n[x,y,z])

    return val


class PolarGridTests(ut.TestCase):

    def test_class(self):

        polar_grid = pg.polarGrid(resolution=res, base_lat_fun=test_fun, RE=RE)

        gridX, gridY, gridZ = polar_grid.get_cartesian_grid()
        gridR, gridL, gridP = polar_grid.get_spherical_grid()

        # Testing for a returned Grid
        self.assertTrue(gridX is not None)
        self.assertTrue(gridY is not None)
        self.assertTrue(gridZ is not None)

        data = np.sqrt(gridX**2 + gridY**2 + gridZ**2)

        d1 = np.isclose(data, RE, atol=1e-12, rtol=0.0)
        d2 = np.isnan(data)
        d3 = d1 | d2

        s1 = np.shape(d1)
        s2 = np.shape(d2)
        s3 = np.shape(d3)

        # Testing self consistency
        self.assertTrue(s1 == s2 == s3, "Failure: {}, {}, {} are not equal".format(s1,s2,s3))
        self.assertTrue(s1 == (1, res[1], res[0]), "Failure: {} != {}".format(s1, [1, res[1], res[0]]))
        self.assertTrue(np.all(d3))

        # Testing lower bounds
        phis = gridP[0,0,:]

        # loop through the phis to check low bound alignment
        for count in range(len(phis)):
            phi = gridP[0,0,count]
            self.assertTrue(phi == phis[count])
            lambdas = gridL[0,:,count]
            fnni = np.argwhere(~np.isnan(lambdas))[0]
            grid_low_lambda = gridL[0,fnni,count]
            grid_phi = gridP[0,fnni,count]

            self.assertTrue(grid_low_lambda == test_fun(phi), msg="Failure in grid lambda alignment")
            self.assertTrue(grid_phi == phi, msg="Failure in grid phi structure Alignment")

        print ("Grid Creation Test Completed Successfully")

    def test_surface_integration(self):
        setup_string2 = setup_string + "\npgrid = pg.polarGrid(resolution={}, " \
                                       "base_lat_fun=lambda longitude: 0)".format(res)

        pgrid = pg.polarGrid(RE=RE, resolution=res, base_lat_fun=test_fun)

        flux = pgrid.get_surface_integral(value_fun=None, value_coord='cart')

        self.assertAlmostEqual(flux, RE**2 * np.pi*2, places=3, msg="Basic hempisphere calculation failed\n"
                                                           "Area calculated: {} -- Area expeceted: {} "
                                                           "-- Difference: {}".format(flux, np.pi*2, RE**2 * np.pi*2 - flux))

        # test flux and L* calculations
        l_tests = np.linspace(start=0, stop=np.pi/2, num=100, endpoint=False)
        for l in l_tests:
            pgrid = pg.polarGrid(resolution=res, base_lat_fun=lambda longitude:l)
            flux = pgrid.get_surface_integral(value_fun=test_value_fun)
            l1 = ic.L_star_from_flux(flux, RE=RE)
            l_ex = ic.analytic_dipole_L(ic.rad_to_deg(l), r=RE)
            self.assertAlmostEqual(l_ex, l1, places=3, msg="L*({}) = {} ==> Expected: {} ==> Diff: {}".format(ic.rad_to_deg(l), l1, l_ex, l_ex-l1))

        print("Surface Integration Test completed Successfully")