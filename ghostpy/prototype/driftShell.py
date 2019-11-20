# File: driftShell.py
# Author: Joshua Murphy
# Project: PhD Dissertation Proposal proof of concept
# ----------------------------------------------------------------------------------
# Purpose:
#       The class in this file is used for finding the drift shell for a given
#       I or K.
# ----------------------------------------------------------------------------------

import bisect as bs
import collections as col
import numpy as np
import paraview.simple as pv
import sys
import vtk

import fieldLine as fl
import inv_common as ih


# TODO: Fix failure to find a drift shell in the realm of the data

class driftShell:
    """
    A class to handle Drift Shells
    """

    def __init__(self, PV_dataObject=None, local_times=None, K=None, b_mirror=None, start_line=None, mode='K', error_tol=None):
        """
        Initialization routine for driftShell class. Sets up the intial configuration of a drift shell
        :param PV_dataObject: Paraview Data object that contains the grid
        :param local_times: array like object of starting local times for finding starting lines (mode='K', K, b_mirror) required.
        :param K: value of K used to define the particle
        :param b_mirror: b_mirror (nT) used to define the particle
        :param start_line: (x,y,z) coordinates for starting line (requires b_mirror)
        :param mode: 'location' or 'K'.
                     'location': starts with a specific location and traces first line there. Uses the given b_mirror.
                     'K': utilizes local time and K to itterate and find the requested b_mirror
        """

        self.data = PV_dataObject
        self.is_valid = True
        self.analytic = False

        if error_tol is None:
            self.error_tol = 1e-5
        else:
            self.error_tol = error_tol

        if mode is 'location_k':
            print "Location (K) based mode"

            if start_line is None or K is None:
                print "ERROR: start_line and K are both required in this mode"
                sys.exit()
            self.start_location = start_line
            self.start_RE = ih.mag(start_line)
            self.K = K
            self.field_lines = col.OrderedDict()
            localtime2 = ih.get_local_time_from_location(start_line)
            self.field_lines[localtime2] = fl.fieldLine(self.data, start=start_line)
            self.b_mirror = self.field_lines[localtime2].get_B_mirror_for_K(self.K)
            for lt in local_times:
                print "Processing Local Time: ", lt
                self.add_field_line_from_K(local_time=lt)
                if not self.field_lines[lt] is None:
                    self.start_RE = ih.mag(self.field_lines[lt].startLoc)

        elif mode is 'location_b':
            print "Location (b_mirror) based mode"

            if start_line is None or b_mirror is None:
                print "ERROR: start_line and b_mirror are both required in this mode"
                sys.exit()
            self.start_location = start_line
            self.start_RE = ih.mag(start_line)
            self.b_mirror = b_mirror
            self.field_lines = col.OrderedDict()
            localtime2 = ih.get_local_time_from_location(start_line)
            self.field_lines[localtime2] = fl.fieldLine(self.data, start=start_line)
            self.K = self.field_lines[localtime2].get_K(self.b_mirror)
            for lt in local_times:
                print "Processing Local Time: ", lt
                self.add_field_line_from_K(local_time=lt)
                if not self.field_lines[lt] is None:
                    self.start_RE = ih.mag(self.field_lines[lt].startLoc)

        elif mode is 'K':
            print "K search based mode"
            if K is None or b_mirror is None or local_times is None:
                print "ERROR: local_times, K, and b_mirror are all required in this mode"
                sys.exit()

            self.start_RE = 5.0
            self.b_mirror = b_mirror
            self.K = K

            self.field_lines = col.OrderedDict()
            for lt in local_times:
                print "Processing Local Time: ", lt
                self.add_field_line_from_K(local_time=lt)
                if not self.field_lines[lt] is None:
                    self.start_RE = ih.mag(self.field_lines[lt].startLoc)

        elif mode is 'analytic_K':
            print "Analytic (K) based mode"
            self.analytic = True
            if start_line is None or K is None:
                print "ERROR: start_line and K are both required in this mode"
                sys.exit()
            self.start_location = start_line
            self.start_RE = ih.mag(start_line)
            self.K = K
            self.field_lines = col.OrderedDict()
            localtime2 = ih.get_local_time_from_location(start_line)
            self.field_lines[localtime2] = fl.fieldLine(start=start_line, error_tol=self.error_tol)
            self.b_mirror = self.field_lines[localtime2].get_B_mirror_for_K(self.K)
            for lt in local_times:
                print "Processing Local Time: ", lt
                self.add_field_line_from_K(local_time=lt)
                if not self.field_lines[lt] is None:
                    self.start_RE = ih.mag(self.field_lines[lt].startLoc)

        else:
            print "Mode ", mode, " is not understood. Exiting."
            sys.exit()

    def add_field_line_from_K(self, local_time):
        """
        Adds a field line object to the drift shell list, based on an iteration over space for K to match B_mirror
        :local_time: the local time for which to search for a field line
        :return: True for success, False for failure
        """
        start_unit = ih.get_location_from_localtime(local_time)
        start_location = list(start_unit * self.start_RE)
        # print start_location
        if not self.analytic:
            newLine = fl.fieldLine(self.data, start=start_location, error_tol=self.error_tol)
        else:
            newLine = fl.fieldLine(start=start_location, error_tol=self.error_tol)
        B_nl = newLine.get_B_mirror_for_K(self.K)
        new_RE = self.start_RE
        low_RE = 1.0
        hi_RE = 3 * new_RE      # TODO: This is a source of problems. needs to fix to outer boundary of grid.

        # print "start_unit: ", start_unit
        # print "start_location: ", start_location
        # print "newLine: ", newLine
        # print "B_nl: ", B_nl
        # print "new_RE: ", new_RE
        # print "low_RE: ", low_RE
        # print "hi_RE:  ", hi_RE
        #
        # print "Start RE: ", self.start_RE

        if B_nl is None:
            print "Line Not Found"
            print "Is Bifurcated? ", newLine.is_bifurcated
            b_values, K_values = ih.dict_to_x_y(newLine.get_k_integrals())
            print "self.K", self.K
            print "Line K Range:   ", K_values[0], ":" ,K_values[-1]
            print "Line B_m Range: ", b_values[0], ":", b_values[-1]

            return False

        while not np.isclose(B_nl, self.b_mirror, atol=self.error_tol, rtol=0.0) and not np.isclose(hi_RE,low_RE, atol=self.error_tol, rtol=0.0):
            dist = np.float64(self.b_mirror-B_nl)
            sep_orig = np.abs((self.b_mirror-B_nl)/self.b_mirror)
            sep = 1.08*sep_orig

            if sep > 1.0 or sep < 0.1:
                print "Fixing SEP: {}".format(sep)
                sep = 0.5

            if B_nl < self.b_mirror or B_nl is None:
                print "Reducing RE"
                hi_RE = new_RE
                new_RE = hi_RE - sep*(hi_RE-low_RE)

            else:
                print "Increasing RE"
                low_RE = new_RE
                new_RE = low_RE + sep*(hi_RE-low_RE)

            print "SEP: {}".format(sep)
            print "DIST: {}".format(dist)
            print "HI - LO: {}-{} = {}".format(hi_RE, low_RE, hi_RE - low_RE)
            print "New RE: {}".format(new_RE)

            new_start = start_unit * new_RE
            if not self.analytic:
                newLine.recompute_field_line(new_start_location=tuple(new_start))
            else:
                newLine.recompute_field_line_analytic(new_start_location=tuple(new_start))

            B_nl = newLine.get_B_mirror_for_K(self.K)
            if B_nl is None:
                print "Re-Adjusting..."

                if new_RE < self.start_RE:
                    while B_nl is None:
                        redux = 0.5 * (hi_RE - new_RE)
                        new_RE += redux
                        new_start = start_unit * new_RE
                        if not self.analytic:
                            newLine.recompute_field_line(new_start_location=tuple(new_start))
                        else:
                            newLine.recompute_field_line_analytic(new_start_location=tuple(new_start))
                        B_nl = newLine.get_B_mirror_for_K(self.K)

                else:
                    while B_nl is None:
                        redux = 0.5 * (new_RE - low_RE)
                        new_RE -= redux
                        new_start = start_unit * new_RE
                        if not self.analytic:
                            newLine.recompute_field_line(new_start_location=tuple(new_start))
                        else:
                            newLine.recompute_field_line_analytic(new_start_location=tuple(new_start))
                        B_nl = newLine.get_B_mirror_for_K(self.K)

                if np.isclose(new_RE, low_RE, atol=self.error_tol, rtol=0.0) or np.isclose(new_RE, self.start_RE, atol=self.error_tol, rtol=0.0):
                    print "ERROR: None Value found"
                    B_nl = None
                    self.is_valid = False
                    break
                else:
                    continue

        if B_nl is None:
            self.field_lines[local_time] = None
            return False
        else:
            self.field_lines[local_time] = newLine
            return True

    def get_phi_locations(self, RE=1):
        phi_locs = dict()

        for key in self.field_lines:
            phi_locs[key] = self.field_lines[key].get_location_for_RE(RE=RE)

        return phi_locs

    def get_new_phi_for_localtime(self, localtime, RE=1):
        """
        returns the intersection location of RE and field line searched for from localtime
        :param localtime:
        :param RE:
        :return:
        """
        # TODO: change this to work with spherical rotation interpolation (great circle)
        start_lines_old = self.get_phi_locations(RE=RE)
        start_lines = dict()
        start_lt = sorted(start_lines.keys())

        # normalize the starting vectors (so we have unit vectors to work with)
        # TODO: Basing this off the keys is wrong, as the field line might be found on a localtime,
        # TODO: but may not intersect anywhere near that line.

        for key in start_lines_old:
            # Need to build new keys based on actual local-time of intersection
            mag_start = start_lines_old[key]/ih.mag(start_lines_old[key])
            new_key = ih.get_local_time_from_location(mag_start)
            start_lines[new_key] = mag_start
            if np.isclose(new_key, 24.0, rtol=0.0, atol=self.error_tol):
                # print "Found 2400: ", new_key
                start_lines[0.0] = mag_start
            elif np.isclose(new_key, 0.0, rtol=0.0, atol=self.error_tol):
                # print "Found 0000: ", new_key
                start_lines[24.0] = mag_start

        # TODO: How to loop if we dont have 2400/0000 as a key?

        start_lt = sorted(start_lines.keys())
        # print start_lt
        index_high = bs.bisect_left(start_lt, localtime)
        index_low = index_high -1

        # adjust to circular
        if index_high is len(start_lt):
            index_high = 0

        latitude_high = np.arcsin(start_lines[start_lt[index_high]][2])
        latitude_low = np.arcsin(start_lines[start_lt[index_low]][2])

        w = (float(localtime) - float(start_lt[index_low])) / (float(start_lt[index_high]) - float(start_lt[index_low]))

        latitude_new = latitude_low + (latitude_high - latitude_low) * w

        loc_localtime = ih.get_location_from_localtime(localtime)

        longitude_new = (np.arctan2(loc_localtime[1],loc_localtime[0]))
        ret = tuple([np.cos(longitude_new) * np.cos(latitude_new)*RE, np.sin(longitude_new) * np.cos(latitude_new)*RE, np.sin(latitude_new)*RE])
        return ret

    def build_polar_cap_grid(self, RE=1, max_lat_divs=18, lon_divs=36, gridType='dict'):
        """
        Generates the grid for integration based on the currently calculated field lines and their interpolations
        :param RE: Distance from center of Earth
        :param max_lat_divs: longitude resolution
        :param lon_divs: latitude resolution
        :param gridType: can be either "dict" or "vtk" or "points"
        :return: eihter a vtkGrid or a dictionary grid based on what is asked for
        """
        grid = col.OrderedDict()
        sph_grid = col.OrderedDict()
        grid_lat = col.OrderedDict()
        grid_lon = col.OrderedDict()

        # latitude spaceing is based on the given latitude
        lon = np.linspace(0, 24, lon_divs, endpoint=False)

        # longitude is based on the lowest latitude.  latitude locations on all other longitude lines use these points
        # let us first get the base latitude/longitudes
        for lt in lon:
            grid[lt]=[]
            sph_grid[lt] = []
            loc = self.get_new_phi_for_localtime(lt, RE)
            # TODO: Fix this hack. We should not be getting invalid values from get new_phi_for_location

            grid_lat[lt] = np.arcsin(loc[2])
            if np.isnan(grid_lat[lt]):
                print "retrieved phi: {}".format(loc)
                print "Error getting Latitude"
                assert False
            grid_lon[lt] = np.arctan2(loc[1], loc[0])

        # find lowest latitude
        low_lat = np.min(np.array(grid_lat.values()))

        # define longitude intervals
        lon_intervals = np.linspace(low_lat, np.pi / 2, max_lat_divs)

        # build the grid:
        pts = vtk.vtkPoints()
        for key in grid_lat:
            for lat in lon_intervals:

                if lat < grid_lat[key]:
                    location = tuple([np.cos(grid_lon[key]) * np.cos(grid_lat[key]) * RE, np.sin(grid_lon[key]) * np.cos(grid_lat[key]) * RE, np.sin(grid_lat[key]) * RE])
                    grid[key].append(location)
                    pts.InsertNextPoint(location)
                    sph_grid[key].append((RE, grid_lon[key], grid_lat[key]))

                else:
                    location = tuple([np.cos(grid_lon[key]) * np.cos(lat)*RE, np.sin(grid_lon[key]) * np.cos(lat)*RE, np.sin(lat)*RE])
                    # We were having problems with duplicates in the integration... hopefully this will help
                    if location not in grid[key]:
                        grid[key].append(location)
                        pts.InsertNextPoint(location)
                        sph_grid[key].append((RE, grid_lon[key], lat))

        # set the points
        vtkGrid = vtk.vtkStructuredGrid()
        vtkGrid.SetPoints(pts)

        if gridType is "dict":
            return grid
        elif gridType is "vtk":
            return vtkGrid
        elif gridType is "points":
            return pts
        elif gridType is "sphereCart":
            return grid, sph_grid
        else:
            print "ERROR: Unknown grid type"
            return None

    def integrate_polar_cap(self, RE=1.0, lat_divs=36, lon_divs=9):
        """
        integrates the polar cap as defined by the current drift trajectory
        :return:
        """
        # print "Integrating Polar Cap"
        grid_cart, grid_sphere = self.build_polar_cap_grid(RE=RE, max_lat_divs=lat_divs, lon_divs=lon_divs, gridType="sphereCart")

        #grid_values = self.get_grid_values(grid_cart)
        grid_values = self.get_dipole_values(grid_cart)

        summation = 0
        summation2 = 0

        # make dpdt pieces
        elements = {}
        elements2 = {}

        for key in grid_sphere:
            elements[key] = []
            elements2[key] = []

        hkeys = grid_sphere.keys()
        num_lon = len(hkeys)

        for index_lon in range(num_lon):
            for index_lat in range(len(grid_sphere[hkeys[index_lon]]) - 1):

                lon = index_lon
                lat = index_lat

                if lon + 1 >= num_lon:
                    lower_lon = num_lon - 1
                    upper_lon = 0
                    # print "Readjusted longitude Key: ", lower_lon, ":", upper_lon

                else:
                    lower_lon = lon
                    upper_lon = lon + 1

                lower_lat = lat
                upper_lat = lat + 1

                # print "lower longitude index: ", lower_lon
                # print "upper longitude index: ", upper_lon
                # print "lower latitude index: ", lower_lat
                # print "upper latitude index: ", upper_lat

                r1 = grid_sphere[hkeys[lower_lon]][lower_lat][0]
                lat2 = grid_sphere[hkeys[lower_lon]][upper_lat][2]
                lat1 = grid_sphere[hkeys[lower_lon]][lower_lat][2]
                lon1 = grid_sphere[hkeys[lower_lon]][lower_lat][1]
                lon2 = grid_sphere[hkeys[upper_lon]][lower_lat][1]

                lat_c = (lat1 + lat2)/2
                lon_c = (lon1 + lon2)/2
                norm_c = ih.get_location_from_theta_phi_r(lat_c, lon_c, RE)

                values1 = grid_values[hkeys[lower_lon]][lower_lat]
                values2 = grid_values[hkeys[lower_lon]][upper_lat]
                values3 = grid_values[hkeys[upper_lon]][lower_lat]
                values4 = grid_values[hkeys[upper_lon]][upper_lat]

                values_c = (np.array(values1) + np.array(values2) + np.array(values3) + np.array(values4))/4

                fluxCenter = np.dot(values_c, norm_c)

                if lon1 < 0:
                    lon1 += 2 * np.pi

                if lon2 < 0:
                    lon2 += 2 * np.pi

                # print "lower longitude, lower latitude (Latitude):", lat1
                # print "lower longitude, upper latitude (Latitude):", lat2
                # print "lower longitude, lower latitude (Longitude):", lon1
                # print "upper longitude, lower latitude (Longitude):", lon2

                # dT needs to be based on co-latitude, not latitude
                dT = lat2 - lat1

                if lon2 < lon1:
                    # print "Updating dP lon2"
                    dP = (lon2 + (2 * np.pi)) - lon1
                else:
                    dP = lon2 - lon1

                area1 = 0.5 * (RE ** 2 * np.cos(lat1) * dT * dP)
                if area1 < 0:
                    if area1 < 0:
                        if np.isclose(0, area1, rtol=0, atol=error_tol):
                            area1 = 0.0
                        else:
                            print "invalid area: {}".format(area1)
                            sys.exit()
                area2 = 0.5 * (RE ** 2 * np.cos(lat2) * dT * dP)
                if area2 < 0:
                    if np.isclose(0, area2, rtol=0, atol=error_tol):
                        area2 = 0.0
                    else:
                        print "invalid area: {}".format(area2)
                        sys.exit()

                elements[hkeys[index_lon]].append((area1+area2) * fluxCenter)



        for key in elements:
            summation += np.sum(elements[key])
            # summation2 += np.sum(elements2[key])

        return summation

    def calculate_L_star(self, RE=1.0, lat_divs=10, lon_divs=10):
        flux = self.integrate_polar_cap(RE=RE, lat_divs=lat_divs, lon_divs=lon_divs)

        return ih.L_star_from_flux(flux=flux, RE=RE)

    def get_grid_values(self, grid):
        gridValues = {}
        for key in grid:
            gridValues[key] = []

        subset = pv.ExtractSubset(Input=self.data)

        # TODO: This is currently hardcoded for the small grid.  Must make it dynamic.
        subset.VOI = [390, 425, 185, 225, 185, 225]

        probe = pv.ProbeLocation(Input=subset, ProbeType='Fixed Radius Point Source')
        probe.ProbeType.Radius = 0.0

        for key in grid:
            for pts in grid[key]:
                probe.ProbeType.Center = pts
                try:
                    value = pv.servermanager.Fetch(probe).GetPointData().GetAbstractArray('B').GetTuple(0)
                    gridValues[key].append(value)
                except:
                    print "Value not found for point: ", pts
                    sys.exit()

        return gridValues

    def get_dipole_values(self, grid):
        dipole_values = {}
        for key in grid:
            dipole_values[key] = ih.get_dipole_value(grid[key])
        return dipole_values
