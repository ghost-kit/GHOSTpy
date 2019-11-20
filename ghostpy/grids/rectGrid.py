# File: rectGrid.py
# Purpose: build a grid of the LFM hexahedral style for populating with dipole or T96 fields
# Project: PhD Thesis
# Date: 30 November 2016

import numpy as np

from t96 import t96_01

from prototype import inv_common as ih


class Grid:
    def __init__(self, extents=None, resolution=None, exclusion_distance=None):
        self._type_ = type

        # check for default dim and extent values
        if resolution is None:
            self._dims_ = [128, 128, 128]
        else:
            self._dims_ = resolution

        if extents is None:
            self._extents_ = [-10, 10, -10, 10, -10, 10]
        else:
            self._extents_ = extents

        if exclusion_distance is None:
            self._exclusion_distance_ = 0.25
        else:
            self._exclusion_distance_ = exclusion_distance

        self.origin_exclusion = []

        # Calculate and populate the grid
        self.x, self.y, self.z = self.build_grid()

    def dipole_field(self, ps=0.):
        """
        Generates a dipole field on the internal grid.
        :param ps: Dipole Tilt. Default is 0.
        :return: X, Y, Z components of the dipole field for the internal grid (GSM Coordinates)
        """
        P = self.x**2
        U = self.z**2
        T = self.y**2
        Q = 30574.0/np.sqrt(P + U + T)**5
        V = 3.0 * self.z * self.x
        sps = np.sin(ps)
        cps = np.cos(ps)

        Bx = Q * ((T+U-2.0*P)*sps-V*cps)
        By = -3.0*self.y*Q*(self.x*sps+self.z*cps)
        Bz = Q * ((P+T-2.0*U)*cps-V*sps)

        for ex in self.origin_exclusion:
            Bx[ex] = np.nan
            By[ex] = np.nan
            Bz[ex] = np.nan

        return Bx, By, Bz

    def t96_field(self, dyn_preassure, dst, imf_by, imf_bz, dipole_tilt):
        """
        Generates the T96 model fields and modifies a dipole
        :param dyn_preassure: Solar Wind Dynamic Pressure
        :param dst: Disturbance (storm time) index.
        :param imf_by: Interplanetary Magnetic Field - Y component (GSM)
        :param imf_bz: Interplanetary Magnetic Field - Z component (GSM)
        :param dipole_tilt:
        :return:
        """
        # get the dipole
        b_x, b_y, b_z = self.dipole_field(ps=dipole_tilt)

        # add the t96 external components
        print "Dst: " , dst
        print "Pdyn: ", dyn_preassure
        print "IMF_By: ", imf_by
        print "IMF_Bz: ", imf_bz
        print "dipole_tilt: ", dipole_tilt

        bx, by, bz = self.t96_wrap(dyn_preassure, dst, imf_by, imf_bz, dipole_tilt)

        b_x = b_x + bx
        b_y = b_y + by
        b_z = b_z + bz

        for ex in self.origin_exclusion:
            b_x[ex] = np.nan
            b_y[ex] = np.nan
            b_z[ex] = np.nan

        return b_x, b_y, b_z

    def build_grid(self):
        # Dimensions

        nx, ny, nz = self._dims_[0], self._dims_[1], self._dims_[2]  # Grid Extent resolution
        lx, ly, lz = self._extents_[0], self._extents_[2], self._extents_[4]
        ux, uy, uz = self._extents_[1], self._extents_[3], self._extents_[5]

        # Coordinates
        X = np.linspace(lx, ux, nx)
        Y = np.linspace(ly, uy, ny)
        Z = np.linspace(lz, uz, nz)

        # Build the Grid
        x = np.zeros((nx, ny, nz))
        y = np.zeros((nx, ny, nz))
        z = np.zeros((nx, ny, nz))

        for k in range(len(Z)):
            for j in range(len(Y)):
                for i in range(len(X)):
                    x[i, j, k] = X[i]
                    y[i, j, k] = Y[j]
                    z[i, j, k] = Z[k]
                    if (-self._exclusion_distance_ < X[i] < self._exclusion_distance_) and (-self._exclusion_distance_ < Y[j] < self._exclusion_distance_) and (-self._exclusion_distance_ < Z[k] < self._exclusion_distance_):
                        if ih.mag([X[i],Y[j],Z[k]]) < self._exclusion_distance_:
                            self.origin_exclusion.append((i, j, k))
                            print "Excluding: ", (i,j,k)

        return x, y, z

    def t96_wrap(self, p_dyn, dst, imf_by, imf_bz, ps):
        parmod = np.array([p_dyn, dst, imf_by, imf_bz, 0, 0, 0, 0, 0, 0])
        print "Calculating the T96 Model... please stand by..."
        t96_v = np.vectorize(t96_01, excluded=['parmod', 'ps'])

        return t96_v(x=self.x, y=self.y, z=self.z, parmod=parmod, ps=ps)