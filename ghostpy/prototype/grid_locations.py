# File: grid_locations.py
# Author: Joshua Murphy
# Project: PhD Dissertation Proposal proof of concept
# ---------------------------------------------------------------------
# Purpose:
#       Lookup table based location finder based on local times.
#       This provides the facilities to have adjustible accuracy based on number of starting points.
#       The system relies on starting local times to generate an interpolated location.
#       All starting locations originate on the equator.
#       provides locations in Cartesian coordinates
# ---------------------------------------------------------------------

import numpy as np

# TODO: class for LUT lookups
#   TODO: self-generating lookup table (dictionary)
#
#   TODO: On setup:   1) Must be provided seed list
#   TODO:             2) initialize LUT with seed list
#
#   TODO: On request: 1) check lookup table.
#   TODO:             2) if not in lookup table, calculate starting point with interpolation
#   TODO:             3) scale starting point to desired RE
#   TODO:             4) Return Location

class lookUpTable:

    def __init__(self, seed_locations=None):
        """
        initialization for the lookUpTable object.  Sets up the lookup table to the seeded local times
        :param seed_locations: local times (floats) in hours (24 hour clock)
        """
        if seed_locations is None:
            seed_locations = [0., 6., 12., 18.]
        self.LUT = {}

        for ang in seed_locations:
            self.LUT[ang] = self._calc_starting_location(ang)
            # print "LT: ", ang, " Loc: ", self.LUT[ang]

    def _calc_starting_location(self, seed):
        """
        calculates a staring location based on circle and local time
        :param seed: Local Time to find unit vector
        :return: unit vector for seeded local time
        """
        clockangle = self._get_clock_angle(seed)
        # need to reverse the x coordinate because... Earth.
        return np.array([-np.cos(clockangle), -np.sin(clockangle), 0.0])

    def _get_clock_angle(self, localtime):
        """
        Generates a clock angle (in radians) from a provided local time (float)
        :param localtime: (float) representing hours on a 24 hour clock (0 - 24 ar valid)
        :return: clock angle in radians
        """
        degs = 0.25 * 60 * localtime
        return degs * np.pi/180
