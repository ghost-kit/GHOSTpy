# File: fieldLine.py
# Author: Joshua Murphy
# Project: PhD Dissertation Proposal proof of concept
# ----------------------------------------------------------------------------------
# Purpose:
#       This file provides field line manipulation functionality for the project.
#       This module requires use of paraview.simple from the paraview application.
# ----------------------------------------------------------------------------------

import bisect as bs
import numpy as np
import paraview.simple as pv
import sys
from collections import OrderedDict
from collections import deque

import inv_common as ih
import valBounds as vb


class fieldLine:
    """ A Class to handle individual field line traces """

    def __init__(self, data, vector_type_name=None, start=None, error_tol=None):
        """
        Constructs the class instance and populates the data channels from Paraview.
        :param data: A Data Object from Paraview (any pipeline data object)
        :param vector_type_name: The name of the VECTOR field you want to trace. Default: ['POINTS', 'B']
        :param start: Starting location for the field-line tracer. Default: [2.0, 0.0, 0.0]
        """
        if vector_type_name is None:
            self.vector = ['POINTS', 'B']
        else:
            self.vector = vector_type_name

        if start is None:
            self.startLoc = (6.6, 0.0, 0.0)
        else:
            self.startLoc = start

        if error_tol is None:
            self.error_tol = 1.0e-5
        else:
            self.error_tol = error_tol

        self.dataObj = data
        self.polar_cap_intersection = dict(North=None, South=None)
        self.is_bifurcated = False
        self.fieldLineData_f = deque()
        self.fieldLineData_b = deque()
        self.fieldLinePoints_f = deque()
        self.fieldLinePoints_b = deque()
        # start the processing
        self.fieldLineObj_f = self.__generate_field_line_obj(direction="FORWARD")
        self.fieldLineObj_b = self.__generate_field_line_obj(direction="BACKWARD")
        self.__extract_field_line_trace()
        self.__resample_field_line()

    #TODO: fix this so that it works in one init
    # def __init__(self, start=None, field_type=None, error_tol=None):
    #     """
    #     create analytic field lines
    #     :param start: Starting location of form [x,y,z]
    #     :param field_type: Valid types are 'dipole', 't96'
    #     """
    #     if field_type is None:
    #         field_type = 'dipole'
    #
    #     if start is None:
    #         self.startLoc = [6.6, 0.0, 0.0]
    #     else:
    #         self.startLoc = start
    #
    #     if error_tol is None:
    #         self.error_tol = 1.0e-5
    #     else:
    #         self.error_tol = error_tol
    #
    #     if field_type is 'dipole':
    #         print "Running in Analytic Dipole Mode"
    #         self.dataObj = None
    #         self.vector = None
    #         self.polar_cap_intersection = dict(North=None, South=None)
    #
    #         # start the processing
    #         path_f, value_f = ih.trace_rk45(x0=self.startLoc, inner_boundary=0.8, h=1e-5, error_tol=self.error_tol, direction="f")
    #         path_b, value_b = ih.trace_rk45(x0=self.startLoc, inner_boundary=0.8, h=1e-5, error_tol=self.error_tol, direction="b")
    #
    #         self.is_bifurcated = False
    #         self.fieldLinePoints_f = deque(path_f)
    #         self.fieldLineData_f = deque(value_f)
    #         self.fieldLinePoints_b = deque(path_b)
    #         self.fieldLineData_b = deque(value_b)
    #
    #         self.__resample_field_line()
    #
    #     elif field_type is 't96':
    #         print "Running in Analytic t96 Mode"
    #         print "Not Yet Implemented"
    #         sys.exit()



    def recompute_field_line(self, new_start_location, error_tol=1e-5):
        """
        Uses the existing field line object and recomputes from a new starting location.
        :param new_start_location: Location on which to re-start the trace
        :return: None
        """
        self.error_tol=error_tol
        self.startLoc = new_start_location
        self.is_bifurcated = False
        self.fieldLineData_f.clear()
        self.fieldLineData_b.clear()
        self.fieldLinePoints_f.clear()
        self.fieldLinePoints_b.clear()
        self.fieldLineObj_f.SeedType.Center = new_start_location
        self.fieldLineObj_b.SeedType.Center = new_start_location
        self.__extract_field_line_trace()
        self.__resample_field_line()

    def recompute_field_line_analytic(self, new_start_location, error_tol=1e-5):
        self.error_tol=error_tol
        self.startLoc = new_start_location
        self.is_bifurcated = False

        path_f, value_f = ih.trace_rk45(x0=self.startLoc, inner_boundary=0.8, h=1e-5, error_tol=self.error_tol, direction="f")
        path_b, value_b = ih.trace_rk45(x0=self.startLoc, inner_boundary=0.8, h=1e-5, error_tol=self.error_tol, direction="b")

        self.is_bifurcated = False
        self.fieldLinePoints_f = deque(path_f)
        self.fieldLineData_f = deque(value_f)
        self.fieldLinePoints_b = deque(path_b)
        self.fieldLineData_b = deque(value_b)

        self.__resample_field_line()

    def __generate_field_line_obj(self, direction='BOTH'):
        """
        Generates the field line object and saves it locally. Requires
        paraview.simple.StreamTracer
        :param direction:
        :return:
        """
        field_line_obj = pv.StreamTracer(Input=self.dataObj, SeedType='Point Source')
        field_line_obj.Vectors = self.vector
        field_line_obj.InterpolatorType = 'Interpolator with Point Locator'
        # self.fieldLineObj.SurfaceStreamlines = 0
        field_line_obj.IntegrationDirection = direction
        field_line_obj.IntegratorType = 'Runge-Kutta 4-5'
        field_line_obj.IntegrationStepUnit = 'Cell Length'
        field_line_obj.InitialStepLength = 0.02
        field_line_obj.MinimumStepLength = 1e-12
        field_line_obj.MaximumStepLength = 0.45
        field_line_obj.MaximumSteps = 100000
        field_line_obj.MaximumStreamlineLength = 60000.0
        field_line_obj.TerminalSpeed = 1e-12
        field_line_obj.MaximumError = self.error_tol
        field_line_obj.ComputeVorticity = 0

        # init the 'Point Source' selected for 'SeedType'
        field_line_obj.SeedType.Center = self.startLoc
        field_line_obj.SeedType.NumberOfPoints = 1
        field_line_obj.SeedType.Radius = 1e-10
        return field_line_obj

    def __extract_field_line_trace(self):
        """
        __extract_field_line_trace(): extracts the field line trace and orders it from center Forward and Backward
                Generates a forward map and a backward map as well as index maps
        :return:    None
        """

        for d in [(self.fieldLineObj_f, self.fieldLineData_f, self.fieldLinePoints_f),
                  (self.fieldLineObj_b, self.fieldLineData_b, self.fieldLinePoints_b)]:
            if d[0] is not None:
                fl_data = pv.servermanager.Fetch(d[0])
                arrays = d[0].PointData.keys()
                if self.vector[1] in arrays:
                    array_points = fl_data.GetPoints()
                    array_data = fl_data.GetPointData().GetAbstractArray(self.vector[1])
                    num_tuples = array_data.GetNumberOfTuples()
                    for i in range(num_tuples):
                        d[1].append(array_data.GetTuple(i))
                        d[2].append(array_points.GetPoint(i))
                else:
                    print "ERROR: ", self.vector[1], " Not found.  Please check script."
                    print "Arrays: ", arrays
                    sys.exit(3)
            else:
                print "D: ", d
                print "ERROR: No data found.  Please check script."
                sys.exit(4)
        return None


    def __resample_field_line(self):
        """
        __resample_field_line(self, direction=0) will re-sample the field line
         in the given direction to match VALUES in the opposite direction.
        :return: None
        """

        mon_inc = lambda L: reduce(lambda a, b: b if a < b else 9999999, L) != 9999999
        mon_dec = lambda L: reduce(lambda a, b: b if a > b else 9999999, L) != 9999999

        # get loc_min
        min_loc, min = self.get_field_min()

        # find the direction minimum exists in
        min_f = min_loc in self.fieldLinePoints_f
        min_b = min_loc in self.fieldLinePoints_b

        min_side_data = None
        min_side_loc = None

        comp_side_data = None
        comp_side_loc = None

        recenter = False

        # Identify group that Bmin resides in
        if min_f and min_b:
            # print "Minimum on BOUNDARY"
            min_side_data = self.fieldLineData_f
            min_side_loc = self.fieldLinePoints_f

            comp_side_data = self.fieldLineData_b
            comp_side_loc = self.fieldLinePoints_b
            recenter = False

        elif min_f and not min_b:
            # print "Minimum in FORWARD group"
            min_side_data = self.fieldLineData_f
            min_side_loc = self.fieldLinePoints_f

            comp_side_data = self.fieldLineData_b
            comp_side_loc = self.fieldLinePoints_b
            recenter = True

        elif min_b and not min_f:
            # print "Minimum in BACKWARD group"
            min_side_data = self.fieldLineData_b
            min_side_loc = self.fieldLinePoints_b

            comp_side_data = self.fieldLineData_f
            comp_side_loc = self.fieldLinePoints_f
            recenter = True

        else:
            print "ERROR - Unknown Problem... No known Bmin"
            sys.exit(5)

        # Recenter trace on Bmin if needed
        # TODO: find faster way of finding index in deque. maybe just walk the list?
        if recenter:
            min_ind = list(min_side_loc).index(min_loc)
            stack_d = []
            stack_l = []
            # Rotate stuff bellow min into stack
            for i in range(min_ind):
                stack_d.append(min_side_data.popleft())
                stack_l.append(min_side_loc.popleft())
            # verify we have rotated correctly
            if min_side_loc[0] is not min_loc:
                print "Error: Not Shifted properly"
                sys.exit(6)

            # remove starting point from compliment side
            comp_side_data.popleft()
            comp_side_loc.popleft()

            # move stack to complement side
            comp_side_data.extendleft(stack_d)
            comp_side_loc.extendleft(stack_l)

            # add minimum to comp side
            comp_side_data.extendleft([min_side_data[0]])
            comp_side_loc.extendleft([min_side_loc[0]])

            # check to make sure we got it correct
            if comp_side_loc[0] is not min_loc and min_side_loc is not min_loc:
                print "ERROR: We did not re-center correctly"
                sys.exit(7)

            # Re-centering complete.  We now need to interpolate the complement side to align values.

        # 2) Truncate traces to 0.5 RE
        min_side_RE = [np.sqrt(i[0]**2 + i[1]**2 + i[2]**2) for i in min_side_loc]
        comp_side_RE = [np.sqrt(i[0]**2 + i[1]**2 + i[2]**2) for i in comp_side_loc]
        min_side_valid_range = [i for i in min_side_RE if i >= 0.5]
        comp_side_valid_range = [i for i in comp_side_RE if i >= 0.5]

        min_to_remove = len(min_side_data) - len(min_side_valid_range)
        cmp_to_remove = len(comp_side_data) - len(comp_side_valid_range)

        for i in range(min_to_remove):
            min_side_data.pop()
            min_side_loc.pop()

        for i in range(cmp_to_remove):
            comp_side_data.pop()
            comp_side_loc.pop()


        # 3) For each Value on min side, find bounding locations for same value on comp side
        loc_bounds = OrderedDict()

        # Center (Bmin)
        Btot_min_side = [np.sqrt(i[0]**2 + i[1]**2 + i[2]**2) for i in min_side_data]
        Btot_comp_side =[np.sqrt(i[0]**2 + i[1]**2 + i[2]**2) for i in comp_side_data]

        self.Btot_min_side = Btot_min_side
        self.Btot_comp_side = Btot_comp_side
        self.Btot_loc_comp_side = [np.sqrt(i[0]**2 + i[1]**2 + i[2]**2) for i in comp_side_loc]
        self.Btot_loc_min_side = [np.sqrt(i[0]**2 + i[1]**2 + i[2]**2) for i in min_side_loc]

        if mon_inc(Btot_comp_side):
            for b in Btot_min_side:
                ls_index = bs.bisect_left(Btot_comp_side, b)
                if ls_index == 0:
                    loc_bounds[b] = [min_side_data[0], vb.valBounds(lowerVal=min_side_data[0],
                                                                    lowerLoc=min_side_loc[0],
                                                                    upperVal=min_side_data[0],
                                                                    upperLoc=min_side_loc[0])]
                else:

                    try:
                        mon_check = [Btot_comp_side[ls_index-1], b, Btot_comp_side[ls_index]]
                        if mon_check[0] <= mon_check[1] <= mon_check[2]:
                            loc_bounds[b] = [min_side_data[Btot_min_side.index(b)],
                                             vb.valBounds(lowerVal=comp_side_data[ls_index-1],
                                                          lowerLoc=comp_side_loc[ls_index-1],
                                                          upperVal=comp_side_data[ls_index],
                                                          upperLoc=comp_side_loc[ls_index])]
                        else:
                            print "Bisection Failure"
                            print mon_check,
                            print mon_inc(mon_check)
                            sys.exit(1)
                    except:

                        if ls_index != len(Btot_comp_side):
                            print "ERROR: Found Bisection Point is in the incorrect location."
                            print "length of Comp Side: ", len(Btot_comp_side)
                            print "bisection index: ", ls_index
                            sys.exit(8)
                        else:
                            # TODO: make sure we are discarding unneeded.
                            # print "Continuing... should I be?"
                            continue

        else:
            print "CAUTION: Field Line B values are non monotonically increasing."



            self.is_bifurcated = True

        # 4) Interpolate across Bounding locations to find location of Value from min side
        # 5) Store new compliment trace

        if not self.is_bifurcated:      # for now, I am just not going to interpolate bifurcated lines.
            b_vals = loc_bounds.keys()
            comp_side_loc.clear()
            comp_side_data.clear()
            for bval in b_vals:
                comp_side_loc.append(loc_bounds[bval][1].get_location(loc_bounds[bval][0]))
                comp_side_data.append(loc_bounds[bval][0])

            # 6) ensure the forward and backward traces are the same length
            longD = None
            longL = None
            shortD = None
            if len(comp_side_data) > len(min_side_data):
                longD = comp_side_data
                longL = comp_side_loc
                shortD = min_side_data
            elif len(min_side_data) > len(comp_side_data):
                longD = min_side_data
                longL = min_side_loc
                shortD = comp_side_data
            else:
                longD = None

            if longD is not None:
                while len(longD) > len(shortD):
                    longD.pop()
                    longL.pop()

            if not set(comp_side_data) == set(min_side_data):
                print "ERROR: Value Lists for Field Line are NOT equivilent... something went wrong..."
                print "Field line that passes through: ", self.startLoc
                for x in range(len(min_side_data)):
                    print  min_side_data[x], " = ", comp_side_data[x]
                print "Length of min side data: ", len(min_side_data)
                print "Length of min side loc:  ", len(min_side_loc)
                print "Length of cmp side data: ", len(comp_side_data)
                print "Length of cmp side loc:  ", len(comp_side_loc)

                sys.exit(9)

        return None

    def get_field_min(self):
        """
        get_field_min() returns the location of the minimum value
        for the field line and its value in the form (loc,val)
        :return: Field Line location and minimum
        """
        # This will eliminate any possibility of loosing the location
        # calculate |B|
        Btot_f = [np.sqrt(component[0]**2 + component[1]**2 + component[2]**2) for component in self.fieldLineData_f]
        Btot_b = [np.sqrt(component[0]**2 + component[1]**2 + component[2]**2) for component in self.fieldLineData_b]

        # get the minimum value and location
        Bmin_f = np.min(Btot_f)
        Bmin_b = np.min(Btot_b)


        if Bmin_f < Bmin_b:
            Bmin_f_index = Btot_f.index(Bmin_f)
            return self.fieldLinePoints_f[Bmin_f_index], Bmin_f
        elif Bmin_b < Bmin_f:
            Bmin_b_index = Btot_b.index(Bmin_b)
            return self.fieldLinePoints_b[Bmin_b_index], Bmin_b
        else:
            Bmin_b_index = Btot_b.index(Bmin_b)
            Bmin_f_index = Btot_f.index(Bmin_f)
            return tuple(self.fieldLinePoints_b[Bmin_b_index]), Bmin_b

    def get_i_integrals(self):
        """
        returns a dictionary containing the b_mirror value as the Key, and I as the value
        :return: dictionary containing {b_mirror: I}
        """
        b_mir_values = self.get_b_mag_values()
        i_dict = {x: self._integrate_i(x) for x in b_mir_values}
        return i_dict

    def get_k_integrals(self):
        """
        returns a dictionary containing the b_mirror value as the key, and K as the value
        :return: dictionary containing {b_mirror: K)
        """
        b_mir_values = self.get_b_mag_values()
        k_dict = {x: self._integrate_i(x) * np.sqrt(x) for x in b_mir_values}
        return k_dict

    def get_b_mag_values(self):
        """
        returns a list of discreete b_mag values in the field trace
        :return: list containing [|B|] contained in trace.
        """
        return [ih.mag(x) for x in self.fieldLineData_f]

    def get_b_mag_forward(self):
        return [ih.mag(x) for x in self.fieldLineData_f]

    def get_b_mag_backward(self):
        return [ih.mag(x) for x in self.fieldLineData_b]

    def _integrate_i(self, B_mirror):
        """
        integrates a specific I for a given B_mirror
        :param B_mirror:
        :return: Value of I
        """
        if not self.is_bifurcated:
            if B_mirror in self.get_b_mag_values():
                # 1) locate B_mirror index
                b_f_mag = [ih.mag(x) for x in self.fieldLineData_f]
                b_mir_ind = b_f_mag.index(B_mirror)

                # 2) Integrate FORWARD to B_mirror
                int_forward = self.trap_int(dataValue=self.fieldLineData_f,
                                            dataLoc=self.fieldLinePoints_f,
                                            start=0, stop=b_mir_ind)

                # 3) Integrate BACKWARD to B_mirror
                int_backward = self.trap_int(dataValue=self.fieldLineData_b,
                                             dataLoc=self.fieldLinePoints_b,
                                             start=0, stop=b_mir_ind)

                # 4) Return BACKWARD + FORWARD
                return int_forward + int_backward

            else:
                b_mags = self.get_b_mag_values()
                hi_ind = bs.bisect_left(b_mags, B_mirror)
                if hi_ind == len(b_mags) or hi_ind == 0:
                    return None
                else:
                    bmir, Iall = ih.dict_to_x_y(self.get_i_integrals())

                    return ih.cub_interp(Iall, bmir, B_mirror)
        else:
            return None

    def has_b_mirror(self, b_mirror):
        """
        returns True if b_mirror is available in this line
        :param b_mirror: the b_mirror being searched for
        :return: True or False
        """
        b_vals = self.get_b_mag_values()
        return b_vals[0] <= b_mirror <= b_vals[-1]

    def get_I(self, b_mirror):
        """
        Returns an interpolated I for a given b_mirror
        :param b_mirror: Value of b_mirror associated with desired I
        :return: Return I if it exists, None otherwise
        """
        return self._integrate_i(b_mirror)

    def has_I(self, I):
        """
        returns true if the field line contains the I needed
        :param I: the I value being checked
        :return: True if I is present in field line, False otherwise
        """
        b_mags = self.get_b_mag_values()
        hi_I = self._integrate_i(b_mags[-1])
        lo_I = self._integrate_i(b_mags[0])

        return lo_I <= I <= hi_I

    def get_K(self, b_mirror):
        """
        gets a K for a given b_mirror on this field line (Interpolated)
        :param b_mirror: value of b_mirror associated with desired K
        :return: returns K if it exists, None otherwise
        """
        I = self.get_I(b_mirror=b_mirror)
        if I is None:
            return None
        else:
            return self.get_I(b_mirror) * np.sqrt(b_mirror)

    def has_K(self, K):
        """
        returns true if K is present in the field line, false otherwise.
        :param K: Value of K being searched for
        :return: True is K is present, False otherwise
        """
        b_mags = self.get_b_mag_values()
        hi_K = self.get_K(b_mags[-1])
        lo_K = self.get_K(b_mags[0])

        return lo_K <= K <= hi_K


    def get_B_mirror_for_I(self, I):
        """
        returns a B_mirror value for a given I.
        :param I: Value of I
        :return: Value of B_mirror associated with I, or None if doesn't exist
        """
        if self.has_I(I):
            b_mirr, I_vals = ih.dict_to_x_y(self.get_i_integrals())
            return ih.cub_interp(b_mirr, I_vals, I)
        else:
            return None

    def get_B_mirror_for_K(self, K):
        """
        returns the B_mirror value for a given K.
        :param K: Value of K
        :return: Value of B_mirror associated with K, or None if doesn't exist
        """
        if self.has_K(K):
            b_mirr, K_vals = ih.dict_to_x_y(self.get_k_integrals())
            return ih.cub_interp(b_mirr, K_vals, K)
        else:
            return None

    def trap_int(self, dataValue, dataLoc, start, stop):
        """
        function utilizes trapezoidal integration to integrate from start Index to stop Index
        :param dataValue: Array Values to integrate
        :param dataLoc: Array of Locations corresponding to Values
        :param start: index to start integration
        :param stop: index to stop integration
        :param f: function to integrate
        :param f_const_args: constant argument list to be passed to function f
        :return: value of 1/2 SUM[k = 1 : N] ( f(k+1) + f(k) ds )
        """
        i_list = []
        for x in range(start, stop):
            try:
                if not self.is_bifurcated:
                    i_list.append((ih.mag(np.array(dataLoc[x + 1]) - np.array(dataLoc[x]))) *
                                  (np.sqrt(1 - ih.mag(dataValue[x+1])/ih.mag(dataValue[stop])) +
                                   np.sqrt(1 - ih.mag(dataValue[x])/ih.mag(dataValue[stop]))))
            except:
                print "ERROR OF SOME KIND"
                print "start: ", start
                print "stop:  ", stop
                print "Current: ", x
                print "Current+1: ", x+1
                print "Size of Data: ", len(dataValue)
                print "Size of Points: ", len(dataLoc)
                # sys.exit(10)

                break

        return 0.5 * np.sum(i_list)

    def get_location_for_B(self, B):

        if self.is_bifurcated:
            print "Bifurcated"
            return None

        b_mags = self.get_b_mag_values()
        b_ind = bs.bisect_left(b_mags, B)

        if b_ind == len(b_mags):
            print "End of List"
            return None

        if b_ind == 0:
            if np.isclose(B, b_mags[0], rtol=0.0, atol=self.error_tol):
                return [self.fieldLinePoints_f[0], self.fieldLinePoints_b[0]]
            else:
                print "b_ind is 0, not same B."
                return None

        bound_f = vb.valBounds(lowerVal=self.fieldLineData_f[b_ind-1],
                               lowerLoc=self.fieldLinePoints_f[b_ind-1],
                               upperVal=self.fieldLineData_f[b_ind],
                               upperLoc=self.fieldLinePoints_f[b_ind])

        bound_b = vb.valBounds(lowerVal=self.fieldLineData_b[b_ind-1],
                               lowerLoc=self.fieldLinePoints_b[b_ind-1],
                               upperVal=self.fieldLineData_b[b_ind],
                               upperLoc=self.fieldLinePoints_b[b_ind])

        loc = [bound_f.get_location(B), bound_b.get_location(B)]
        return loc

    def get_location_for_RE(self, RE):
        """
        Return the coordinates on the Line where it crosses the given RE distance from Earth.  Returns None if
        no such crossing occurs.
        :param RE: Distance (in RE) for where we want the coordinates on the line.
        :return: list of x,y,z coordinates for where the line crosses the RE distance specified. None otherwise.
        """

        # Get the line in RE format
        RE_f = [ih.mag(x) for x in self.fieldLinePoints_f]

        #if was registerd as bifrucated, we have to do a little extra work
        if self.is_bifurcated:
            print "Caution: This line is registering as bifurcated.\nAttempting to find the intersection point"
            print "monotonically decreasing test: {}".format(ih.mon_dec(RE_f))
            if not ih.mon_dec(RE_f):
                print "Need to adjust trace"
                print "exporting Data Dump"
                pos = np.asarray(self.fieldLinePoints_f)
                dat = np.asarray(self.fieldLineData_f)
                dpos = np.diff(self.fieldLinePoints_f, axis=0)
                ddat = np.diff(self.fieldLineData_f, axis=0)

                np.savetxt("errors/pos_{}.csv".format(self.startLoc), pos, delimiter=",")
                np.savetxt("errors/dat_{}.csv".format(self.startLoc), dat, delimiter=",")
                np.savetxt("errors/dpos_{}.csv".format(self.startLoc), dpos, delimiter=",")
                np.savetxt("errors/ddat_{}.csv".format(self.startLoc), dpos, delimiter=",")



        ind_f = None
        loc_f = None

        # reverse the list
        RE_f.reverse()

        if RE < RE_f[0] or RE > RE_f[-1]:
            print "Out of bounds:"
            print "RE: ", RE
            print "RE_f: ", RE_f[0], ":", RE_f[-1]
        else:
            ind_f = bs.bisect_left(RE_f, RE)
            # print "Index: ", ind_f

            bound_f = vb.valBounds(lowerVal=RE_f[ind_f-1],
                                   lowerLoc=self.fieldLinePoints_f[-ind_f],
                                   upperVal=RE_f[ind_f],
                                   upperLoc=self.fieldLinePoints_f[-(ind_f+1)])

            wa=0.0
            loc_f = bound_f.get_location(RE, weight_adjust=wa)
            loc_RE = ih.mag(loc_f)

            # TODO: make this a percentage... or don't do it at all, as this can cause an infinite loop like behaviour.
            # while not np.isclose(loc_RE, RE):
            #     if loc_RE < RE:
            #         wa += 0.0001
            #         loc_f = bound_f.get_location(RE, weight_adjust=wa)
            #         loc_RE = ih.mag(loc_f)
            #         print "Adjusting: ", wa,
            #         print "new RE: ", loc_RE
            #
            #     if loc_RE > RE:
            #         wa -= 0.0001
            #         loc_f = bound_f.get_location(RE, weight_adjust=wa)
            #         loc_RE = ih.mag(loc_f)
            #         print "Adjusting: ", wa,
            #         print "new RE: ", loc_RE

        # print "Location for ", RE, " RE: ",loc_f
        # print "Verification: ", RE, " = ", ih.mag(loc_f)

        # TODO: calculate and return the backward location as well. (if needed)

        return loc_f





