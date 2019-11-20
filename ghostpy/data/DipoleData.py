# Creates a GHOSTpy data object for dipole field

import numpy as np

import GpData as gpd
from ghostpy.algorithms import DipoleField as dpf
from ghostpy.algorithms import common as algc
from ghostpy.algorithms import FieldTracers as fts

class DipoleData(gpd.data):
    """
    The object provides dipole data when instantiated.   The default will provide a dipole for Earth with no dipole tilt.
    """

    def __init__(self, tilt=0, Be=3.15e4, name="dipole"):
        self.tilt = tilt
        self.Be = Be
        self.name = name
        self.inner_boundary = 1.0
        self.trace_boundary = 0.8
        self.integrator = fts.SPrk45(self)

    def __str__(self):
        return "GHOSTpy {} data object".format(self.name)

    def get_name(self):
        return self.name

    def get_xyz(self, xyz):
        x = xyz[0]
        y = xyz[1]
        z = xyz[2]

        return self.__dipole_values__(x, y, z)

    def get_rlp(self, rlp):
        xyz = algc.sphere_to_cart(r=rlp[0], lam=rlp[1], phi=rlp[2])
        return self.get_xyz(xyz)

    def __dipole_values__(self, x, y, z):
        return np.array(dpf.dipole_field(x=x, y=y, z=z, ps=self.tilt, Be=self.Be))

    def set_inner_boundary(self, re=1.0):
        self.inner_boundary = re

    def get_calc_boundary(self):
        return self.inner_boundary

    def get_trace_boundary(self):
        return self.trace_boundary

    def get_integrator(self):
        return self.integrator
