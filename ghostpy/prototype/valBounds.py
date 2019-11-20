# File: valBounds.py
# Author: Joshua Murphy
# Project: PhD Dissertation Proof of Concept
# -----------------------------------------------------------------------------------
# Purpose:
#       Class to aid in resampling. Holds bounds and their values, calculates
#       location of resample between the two based on linear interp.
#
# -----------------------------------------------------------------------------------
from inv_common import *
np.seterr(invalid='ignore')


class valBounds:

    def __init__(self, lowerVal, lowerLoc, upperVal, upperLoc):
        self.lowerBoundVal = np.array(lowerVal)
        self.lowerBoundLoc = np.array(lowerLoc)
        self.upperBoundVal = np.array(upperVal)
        self.upperBoundLoc = np.array(upperLoc)
        self.w = None

    def get_weight(self, newVal):
        if mag(self.lowerBoundVal) <= mag(newVal) <= mag(self.upperBoundVal):
            self.w = (mag(newVal) - mag(self.lowerBoundVal)) / (mag(self.upperBoundVal) - mag(self.lowerBoundVal))
            if np.isnan(self.w):
                self.w = 0
        else:
            self.w = None
            return None

    def get_b_vector(self, newVal):
        if self.w is None:
            self.get_weight(newVal=newVal)
            if self.w is None:
                return None

        d = (self.upperBoundVal - self.lowerBoundVal) * self.w
        b_vec = tuple(self.lowerBoundVal + d)

        # print("w={}\nd={}\nlowerBoundVal={}".format(self.w,d, self.lowerBoundVal))
        # print("newVal: {}\nIntVal: {}".format(mag(newVal), mag(b_vec)))

        # This assertion is to check that we are within a tenth of a nano-tesla of the desired quantity.  we may be able to do better somehow... need to think about the interpolation some
        # print ("b_vec: {}".format(b_vec))
        assert np.isclose(mag(b_vec), mag(newVal), rtol=0.0, atol=1.5), "Getting b vector failed... newVal: {} ... magnitude of retrieved vector: {}".format(mag(newVal), mag(b_vec))
        return b_vec

    def get_location(self, newVal):

        if self.w is None:
            self.get_weight(newVal=newVal)
            if self.w is None:
                return None

        d = (self.upperBoundLoc - self.lowerBoundLoc) * self.w
        xyz = tuple(self.lowerBoundLoc + d)
        return xyz
