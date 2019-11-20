import numpy as np
import scipy.interpolate as sint

import ghostpy.algorithms.FieldTracers as ft
import ghostpy.algorithms.common as algc
import ghostpy.data.GpData as gpd
from scipy.signal import savgol_filter



# TODO NOTES: 27 Jan 2017: Utilize Spline Interpolation from scipy.interpolate.spline for K and B calcs.

class FieldLine:
    """
    FieldLIne: this class is a container to hold information regarding the 2nd adiabatic invariant of geomagnetically trapped particle motion.
    """

    def __init__(self, data=None, start=None, error_tol=None, smooth=19):
        """
        Intializes the field line object

        :param data: A valid data source is a ghostpy data object, and is required to provide data to the field line tracer
            :type data: gpd.data
        :param start:  This is the staring location for the field line trace in "X,Y,Z" cartesian coordinates. (RE units)
            :type start: np.array
        :param error_tol:   The Error Tolerance for the field line tracing. This will default to 1E-5
            :type error_tol: float

        """

        # initialize internal variables
        self.smooth = smooth
        self.start = None
        self.start_idx = None
        self.data = None
        self.error_tol = None
        self.trace_n = None
        self.trace_s = None
        self.m_trace = None
        self.m_trace_re = None
        self.m_trace_data = None
        self.m_trace_b_mirror = None
        self.trace_data_n = None
        self.trace_data_s = None
        self.aligned = False
        self.bifurcated = False
        self.resampled = False
        self.valid = True  # this is not the best way to do this, but I need to get it working...
        self.valid_code = 0
        self.integrator = data.get_integrator()
        self.K = None
        self.sK = None  # sign of diff(K)
        self.sKr = None # rolled to get previous value
        self.B = None
        self.sB = None  # sign of diff(B)
        self.sBr = None # Rolled to get previous value
        # TODO: check assertions for proper algorithm arguments

        assert isinstance(data, gpd.data)
        self.data = data
        self.safe_boundry = 3.1 #self.data.get_calc_boundary()*1.1


        if error_tol is None:
            self.error_tol = 1e-6
        else:
            self.error_tol = error_tol

        assert self.integrator is not None, "An integrator must be specified"

        self.set_start_location(start=start)

    def __str__(self):
        return "GHOSTpy FieldLine Object"

    def __repr__(self):
        return "GHOSTpy: Field Line Tracing Object\n\tdata: {}\n\tstart: {}\n" \
               "\ttrace_mode: {}\n\terror_tol: {}".format(self.data.get_name(), self.start, self.data.get_name(),
                                                          self.error_tol)

    def clear(self):
        self.start = None
        self.start_idx = None
        self.trace_n = None
        self.trace_s = None
        self.m_trace = None
        self.m_trace_re = None
        self.m_trace_data = None
        self.m_trace_b_mirror = None
        self.trace_data_n = None
        self.trace_data_s = None
        self.aligned = False
        self.bifurcated = False
        self.resampled = False
        self.valid = True  # this is not the best way to do this, but I need to get it working...
        self.valid_code = 0
        self.K = None
        self.sK = None  # sign of diff(K)
        self.sKr = None # rolled to get previous value
        self.B = None
        self.sB = None  # sign of diff(B)
        self.sBr = None # Rolled to get previous value



    def get_trace_mode(self):
        """
        returns the name of the Data object being used in the trace
        :return: string representing the trace mode, as is defined by the data-reading plugin
        """
        assert self.data is not None
        return self.data.get_name()

    def set_start_location(self, start):
        """
        sets the starting location for the trace
        :param start: sets a new starting location and calculates the new line
        :return: None
        """
        self.clear()
        assert start is not None, "You must specify a starting location for the field line trace"
        assert len(start) == 3
        self.start = np.array(start)
        self.__trace_line__()
        self.__kb_model__()

    def get_start_location(self):
        """
        returns the currently defined starting location of the field line trace
        :return: Starting location vector (x,y,z)
        """
        return self.start

    def get_footprint_xyz(self, re=1.0):
        if self.trace_n is None:
            if self.trace_s is not None:
                if algc.mag(self.trace_s[0]) < self.safe_boundry:
                    return self.trace_s[0]
            return [np.NaN, np.NaN, np.NaN]
        assert self.trace_n is not None, "ERROR: Cannot Find North Trace Data."
        rd = algc.mag(self.trace_n)

        # print ("Radial Distance: {}".format(rd))
        # print ("North Trace: {}".format(self.trace_n))
        # print ("Looking for RE: {}".format(re))
        assert self.trace_n is not None, "ERROR: Cannot Find North Trace Data."
        rd = algc.mag(self.trace_n)
        # new_pos = sint.spline(xk=rd, yk=self.trace_n, xnew=re, order=1)
        if not np.any((rd < re)):
            re = np.min(rd)
        new_pos = algc.lin_interp(v=re, X=rd, F=self.trace_n)
        # print ("Footprint from FL: {}".format(algc.cart_to_sphere(new_pos)))
        # print ("Radial Distance: {}".format(rd))
        # print ("Monotonic? {}".format(algc.mon_dec(rd) or algc.mon_inc(rd)))
        # print ("RE requested: {}".format(re))
        if new_pos is None:
            if self.m_trace is not None and self.m_trace_re[0] < self.safe_boundry:
                return self.m_trace[0]
            return [np.NaN, np.NaN, np.NaN]
        return new_pos.flatten()

    def get_footprint_rlp(self, re=1.0):
        if self.trace_n is None:
            if self.trace_s is not None:
                if algc.mag(self.trace_s[0]) < self.safe_boundry:
                    return self.trace_s[0]
            return [np.NaN, np.NaN, np.NaN]
        assert self.trace_data_n is not None, "ERROR: Cannot Find North Trace Data."
        xyzf = self.get_footprint_xyz(re=re)
        rlpf = np.array(algc.cart_to_sphere(xyzf))
        # print ("Footprint RLPF: {}".format(rlpf))
        return rlpf

    def get_b(self, k=0):
        b = self.__get_b_mirror__(k=k)
        # print("B: {}".format(b))
        return b


    def __trace_line__(self):
        """
        internal method to perform the line trace
        :return: None -- internal method
        """
        self.trace_n, self.trace_data_n = self.integrator.integrate(x0=self.start, direct='f', error_tol=self.error_tol)
        self.trace_s, self.trace_data_s = self.integrator.integrate(x0=self.start, direct='b', error_tol=self.error_tol)

        if self.trace_n is None or self.trace_s is None:
            self.valid = False
            self.valid_code = -1
            return

        re_n = algc.mag(self.trace_n)
        re_s = algc.mag(self.trace_s)


        if len(self.trace_n) > 0:
            re_n_1 = re_n[0]
            re_n_2 = re_n[-1]
            re_n_max = np.nanmax(re_n)
            if len(re_s) > 0:
                re_s_max = np.nanmax(re_s)
            else:
                re_s_max = np.nan

        else:
            re_n_1 = 1000.
            re_n_2 = 1000.
            re_n_max = np.nan
            re_s_max = np.nan

        if len(self.trace_s) > 0:
            re_s_1 = algc.mag(self.trace_s[0])
            re_s_2 = algc.mag(self.trace_s[-1])
            re_s_max = np.nanmax(re_n)
            if len(re_n) > 0:
                re_n_max = np.nanmax(re_n)
            else:
                re_n_max = np.nan

        else:
            re_s_1 = 1000.
            re_s_2 = 1000.
            re_n_max = np.nan
            re_s_max = np.nan


        if re_n_1 <= self.safe_boundry and re_n_2 <= self.safe_boundry and re_n_max > re_s_max:
            # Full trace in North... must flip
            print ("Full Trace North: {}".format(len(self.trace_n)))
            self.start_idx = len(self.trace_n)-1
            self.m_trace = np.flipud(self.trace_n)
            self.m_trace_data = np.flipud(self.trace_data_n)
            self.m_trace_b_mirror = algc.mag(self.m_trace_data)
            self.m_trace_re = algc.mag(self.m_trace)

            print ("Full Trace RE:\n{}".format(self.m_trace_re))

        elif re_s_1 <= self.safe_boundry and re_s_2 <= self.safe_boundry and re_s_max > re_n_max:
            print ("Full Trace South")
            # Full trace in South... no flip needed
            self.start_idx = 0
            self.m_trace = self.trace_s
            self.m_trace_data = self.trace_data_s
            self.m_trace_b_mirror = algc.mag(self.m_trace_data)
            self.m_trace_re = algc.mag(self.m_trace)

        elif re_n_2 <= self.safe_boundry and re_s_2 <= self.safe_boundry:
            # print ("Combined Trace")
            # Full trace in combination... must combine
            self.start_idx = len(self.trace_n) - 1

            data_array = np.delete(self.trace_data_s, 0, axis=0)
            data_array = np.concatenate([np.flipud(self.trace_data_n), data_array], axis=0)

            # Combine North and South Location Arrays
            # Values should move from north to south along the line
            loc_array = np.delete(self.trace_s, 0, axis=0)
            loc_array = np.concatenate([np.flipud(self.trace_n), loc_array], axis=0)

            self.m_trace = loc_array
            self.m_trace_data = data_array
            self.m_trace_b_mirror = algc.mag(data_array)
            self.m_trace_re = algc.mag(loc_array)

        else:
            self.valid = False
            self.valid_code = -2
            return

        if self.smooth > 0:
            try:
                # print ("heavy")
                self.m_trace_b_mirror = savgol_filter(self.m_trace_b_mirror, self.smooth, 2)
            except TypeError:
                pass

    def __resample_line__(self, debug=False):
        # TODO: Remove all clauses for resample

        self.resampled = True

    def __get_master_trace__(self):

        return self.m_trace, self.m_trace_data


    def get_b_loc(self, b):
        if self.B is None:
            return None

        ba = self.B
        bb = np.roll(ba, 1)
        bb[0] = np.NaN
        bc = np.roll(ba, -1)
        bc[-1] = np.NaN

        pts1a = ba <= b
        pts1b = b < bc
        pts1 = np.logical_and(pts1a, pts1b)

        # pts2a = bb < b
        # pts2b = b <= ba
        # pts2 = np.logical_and(pts2a, pts2b)

        pts3a = bc < b
        pts3b = b <= ba
        pts3 = np.logical_and(pts3a, pts3b)

        # pts4a = ba <= b
        # pts4b = b < bb
        # pts4 = np.logical_and(pts4a, pts4b)

        # pts = np.logical_or(pts1, pts2)
        pts = np.logical_or(pts1, pts3)
        # pts = np.logical_or(pts, pts4)

        wpts = np.argwhere(pts).flatten()

        return wpts

    def get_k(self, b_mirror):
        """
        internal method to calculate K for the given b_mirror
        :param b_mirror: float: b_mirror value for calculating K
        :return: value of K for given b_mirror - internal method
        """

        klist = self.__get_k_list__(b_mirror)
        print ("K-list: {}".format(klist))

        if len(klist) > 0:
            return np.max(klist)
        else:
            return np.NaN

        # # print ("Valid Line? {} = {}".format(self.valid, self.valid_code))
        #
        # if self.B is None:
        #     return -1
        #
        # idx1 = (np.abs(self.B - b_mirror)).argmin()
        #
        # if np.isclose(self.B[idx1], b_mirror, rtol=0, atol=1e-12):
        #     print ("Found VERY close (or exact) value for K")
        #     return self.K[idx1]
        #
        # sb = self.sK < 0  # Signals if next value is lower than current one
        # sbr = self.sKr > 0  # Signals if the previous value is lower than current one
        #
        # sb2 = self.sK > 0  # Signnals if next value is increasing.
        # sbr2 = self.sK < 0  # Signals if previous value is larger than current value
        #
        # fpoint = self.B[idx1]
        #
        # if fpoint > b_mirror:
        #     # print("use sb/sbr")
        #     if idx1 >= len(sb) or (not sb[idx1] and not sbr[idx1]):
        #         print ("Stuck in local minimum or not available on line")
        #         return -1
        #
        #     # print ("sb[idx]".format(sb[idx1]))
        #     if sb[idx1]:
        #         # print("Can interpolate ahead")
        #         lidx = idx1 + 1
        #         hidx = idx1
        #     if sbr[idx1]:
        #         # print("Can interpolate behind")
        #         lidx = idx1 - 1
        #         hidx = idx1
        # else:
        #     # print("use sb2/sbr2")
        #     if idx1 >= len(sb2) or (not sb2[idx1] and not sbr2[idx1]):
        #         # print ("Stuck in local minimum or not available on line")
        #         return -1
        #
        #     if sb2[idx1]:
        #         # print("Can interpolate ahead")
        #         hidx = idx1 + 1
        #         lidx = idx1
        #     if sbr2[idx1]:
        #         # print("Can interpolate behind")
        #         hidx = idx1 - 1
        #         lidx = idx1
        #
        # hb = self.B[hidx]
        # lb = self.B[lidx]
        #
        # assert hb > b_mirror > lb, "Bounding Issue when finding K"
        #
        # # print ("HB: {}".format(hb))
        # # print ("RB: {}".format(b_mirror))
        # # print ("LB: {}".format(lb))
        #
        # # interpolation weight
        # lw = (b_mirror - lb) / (hb - lb)
        # # print ("Interpolation Weight: {}".format(lw))
        #
        # k = self.K[lidx] + (self.K[hidx] - self.K[lidx]) * lw
        # # print ("K: {}".format(k))
        #
        # return k

    def __get_k_list__(self, b):
        if self.B is None:
            return None

        ka = self.K
        kb = np.roll(ka, 1)
        kb[0] = np.NaN
        kc = np.roll(ka, -1)
        kc[-1] = np.NaN

        ba = self.B
        bb = np.roll(ba, 1)
        bb[0] = np.NaN
        bc = np.roll(ba, -1)
        bc[-1] = np.NaN

        pts1a = ba <= b
        pts1b = b < bc
        pts1 = np.logical_and(pts1a, pts1b)

        pts2a = bb < b
        pts2b = b <= ba
        pts2 = np.logical_and(pts2a, pts2b)

        pts3a = bc < b
        pts3b = b <= ba
        pts3 = np.logical_and(pts3a, pts3b)

        pts4a = ba <= b
        pts4b = b < bb
        pts4 = np.logical_and(pts4a, pts4b)

        pts = np.logical_or(pts1, pts2)
        pts = np.logical_or(pts, pts3)
        pts = np.logical_or(pts, pts4)

        wpts = np.argwhere(pts).flatten()

        bav = ba[wpts]
        kav = ka[wpts]

        bbv = bb[wpts]
        kbv = kb[wpts]

        bcv = bb[wpts]
        kcv = kc[wpts]

        klist = []

        for idx in range(len(wpts)):
            bv = np.array([bav[idx], bbv[idx], bcv[idx]])
            kv = np.array([kav[idx], kbv[idx], kcv[idx]])

            # print ("idx: {}".format(idx))
            li = np.argmin(bv)
            # print (li)
            ui = np.argmax(bv)
            # print (ui)

            lbv = bv[li]
            ubv = bv[ui]
            lkv = kv[li]
            ukv = kv[ui]

            bstat = b <= bav[idx]
            # print (bstat)

            if bstat:
                lbv = bv[li]
                lkv = kv[li]
                ubv = bav[idx]
                ukv = kav[idx]

            else:
                ubv = bv[ui]
                ukv = kv[ui]
                lbv = bav[idx]
                lkv = kav[idx]

            # print ("LKV: {}".format(lbv))
            # print ("K: {}".format(k))
            # print ("UKV: {}".format(ubv))

            if ubv - lbv == 0:
                lw = 0.0
            else:
                lw = (b - lbv) / (ubv - lbv)

            if np.abs(lw) <= 1.0:
                k = kv[li] + (kv[ui] - kv[li]) * lw
                klist.append(k)
            else:
                pass
        klist = np.array(sorted(klist))

        return klist

    def __get_b_mirror_list__(self, k):
        if self.K is None:
            return None

        ba = self.B
        bb = np.roll(ba, 1)
        bb[0] = np.NaN
        bc = np.roll(ba, -1)
        bc[-1] = np.NaN

        ka = self.K
        kb = np.roll(ka, 1)
        kb[0] = np.NaN
        kc = np.roll(ka, -1)
        kc[-1] = np.NaN

        pts1a = ka <= k
        pts1b = k < kc
        pts1 = np.logical_and(pts1a, pts1b)

        pts2a = kb < k
        pts2b = k <= ka
        pts2 = np.logical_and(pts2a, pts2b)

        pts3a = kc < k
        pts3b = k <= ka
        pts3 = np.logical_and(pts3a, pts3b)

        pts4a = ka <= k
        pts4b = k < kb
        pts4 = np.logical_and(pts4a, pts4b)

        pts = np.logical_or(pts1, pts2)
        pts = np.logical_or(pts, pts3)
        pts = np.logical_or(pts, pts4)

        wpts = np.argwhere(pts).flatten()

        kav = ka[wpts]
        bav = ba[wpts]

        kbv = kb[wpts]
        bbv = bb[wpts]

        kcv = kb[wpts]
        bcv = bc[wpts]

        blist = []

        for idx in range(len(wpts)):
            kv = np.array([kav[idx], kbv[idx], kcv[idx]])
            bv = np.array([bav[idx], bbv[idx], bcv[idx]])

            # print ("idx: {}".format(idx))
            li = np.argmin(kv)
            # print (li)
            ui = np.argmax(kv)
            # print (ui)

            lkv = kv[li]
            ukv = kv[ui]
            lbv = bv[li]
            ubv = bv[ui]

            kstat = k <= kav[idx]
            # print (kstat)

            if kstat:
                lkv = kv[li]
                lbv = bv[li]
                ukv = kav[idx]
                ubv = bav[idx]

            else:
                ukv = kv[ui]
                ubv = bv[ui]
                lkv = kav[idx]
                lbv = bav[idx]

            # print ("LKV: {}".format(lkv))
            # print ("K: {}".format(k))
            # print ("UKV: {}".format(ukv))

            if ukv-lkv == 0:
                lw = 0.0
            else:
                lw = (k-lkv)/(ukv-lkv)

            if np.abs(lw) <= 1.0:
                # print ("LW: {}".format(lw))
                b = bv[li] + (bv[ui] - bv[li]) * lw
                # print ("B: {}".format(b))

                blist.append(b)
            else:
                # print("KV: {}".format(kv))
                # print("BV: {}".format(bv))
                pass
        blist = np.array(sorted(blist))
        # print ("B List: {}".format(blist))

        return blist

    def __get_min_b_gap__(self, k, b):
        gap_min = None
        b_list = self.__get_b_mirror_list__(k=k)

        if b_list is not None:
            gap_list = b_list - b
            try:
                gap_min = np.nanmin(gap_list)
                # gap_min = np.nanmax(gap_list)
                # gap_min = np.nanmedian(gap_list)
                # gap_min = np.mean(gap_list)

            except ValueError:
                pass

        return gap_min

    def __get_b_mirror__(self, k):

        """
        calculates a b_mirror for a given k.
        :param k: The K value for which we are looking to find an associated b_mirror. (A B_mirror for the given K may not exist)
        :return: Return b_mirror for the given K on the field line, or None if it does not exist.
        """
        try:
            b = np.max(self.__get_b_mirror_list__(k=k))
        except ValueError:
            return -1

        if b is None:
            b = -1

        return b

    def __min_B__(self):
        """
        returns the minimum values for the northern and southern hemispheres
        :return: 2 x np.array: min_n, min_s of form [index value].  Should cast index to int prior to use
        """
        assert self.trace_n is not None and self.trace_data_s is not None, "ERROR: Field Line Trace is of None value. Starting Location: {}".format(
            self.start)
        assert self.trace_data_n is not None and self.trace_data_s is not None, "ERROR: No field line trace data. Starting Location: {}".format(
            self.start)

        if len(self.trace_n) == 0 or len(self.trace_data_n) == 0:
            # print ("Invalid Trace")
            # print ("Zero Length Trace Produced.")
            # print ("North Trace: {}".format(self.trace_n))
            # print ("South Trace: {}".format(self.trace_s))
            return None, None

        b_mag_n = algc.mag(self.trace_data_n)
        b_mag_s = algc.mag(self.trace_data_s)

        min_n = np.min(b_mag_n)
        min_n_i = np.where(b_mag_n == min_n)[0][0]
        min_s = np.min(b_mag_s)
        min_s_i = np.where(b_mag_s == min_s)[0][0]

        return np.array([min_n_i, min_n]), np.array([min_s_i, min_s])

    def __all_k__(self, debug=False):
        """
        calculates all K for the b_mirror values in the trace
        :return: True if values generated, False otherwise
        """
        if self.K is None:
            return self.__kb_model__()

    def __kb_model__(self):

        if self.m_trace_b_mirror is None:
            return False

        n_k = []
        n_b = []

        td = self.m_trace_b_mirror.copy()
        tl = self.m_trace.copy()


        # for idx in range(len(self.m_trace_b_mirror)):
        bf, kf = self.__k_mod__(mdata=td, mtrace=tl)
        bb, kb = self.__k_mod__(mdata=np.flipud(td), mtrace=np.flipud(tl))
        kb = np.flipud(kb)
        n_k = np.nanmax([kf,kb], axis=0)

        n_b = bf

        n_k[[0,-1]] = np.NaN


        self.K = np.array(n_k)
        self.B = np.array(n_b)

        # print (self.K)

        #
        # n_k = []
        # n_b = []
        # for idx in range(len(self.m_trace_b_mirror)):
        #     nk, nb = self.__get_k_b_for_fl_trace_pt__(idx=idx)
        #     # print ("IDX: {}".format(idx))
        #     # print ("NK: {}".format(nk))
        #     # print ("NB: {}".format(nb))
        #
        #     if nk != -1:
        #         n_k.append(nk)
        #     else:
        #         n_k.append(np.NaN)
        #
        #     n_b.append(nb)
        #
        # self.K = np.array(n_k, dtype=np.float64)
        # self.B = np.array(n_b, dtype=np.float64)
        #
        # if self.K[0] == 0.0:
        #     self.K[0] = np.NaN
        #
        # if self.K[-1] == 0.0:
        #     self.K[-1] = np.NaN

        return True


    def __k_mod__(self, mdata, mtrace):

        assert mdata is not None and mtrace is not None

        bdata = mdata.copy()
        bdata2 = bdata.copy()
        mtrace2 = mtrace.copy()
        re_north = algc.mag(mtrace2[0])
        re_south = algc.mag(mtrace2[-1])
        bm_north = bdata[0]
        bm_south = bdata[-1]
        calc_bound = self.data.get_calc_boundary() * 1.5
        len_data = len(mtrace2)
        if (not np.isclose(re_north, calc_bound) and re_north > calc_bound) or (not np.isclose(re_south, calc_bound) and re_south > calc_bound ):
            # print ("Field line through starting point {} does not return to boundary.".format(self.start))
            # print ("Required Boundary: {} RE".format(calc_bound))
            # print ("Line Bounds: Start: {} RE, End: {} RE".format(re_north, re_south))
            self.valid_code = -2
            self.valid = False
            return None, None

        # list for K:
        klist = []
        blist = []

        start = 0
        while bm_north > bm_south:
            bm_north = bdata[start]
            klist.append(np.NaN)
            blist.append(bdata2[start])
            start += 1

        for idx in np.arange(start=start, stop=len_data):
            cb = bdata2[idx]

            # find the bounding point index.  If is same as idx, then no bounding point exists
            try:
                bidx = np.argmax(bdata[idx + 1:] > cb) + idx + 1
            except ValueError:
                klist.append(np.NaN)
                blist.append(bdata2[idx])
                break

            int_list = bdata[idx:(bidx+1)]
            loc_array = mtrace2[idx:(bidx + 1)]

            # interpolate the last value
            bm = cb
            den = (int_list[-1] - int_list[-2])

            if den == 0:
                print ("Points contain equivilant b_mirrors")
                print ("B_mirror: {}".format(bm))
                print ("Setting Weight to 0")

                lw = 0
            else:
                lw = (bm - int_list[-2]) / den

            bend = int_list[-2] + ((int_list[-1] - int_list[-2]) * lw)
            lend = loc_array[-2] + ((loc_array[-1] - loc_array[-2]) * lw)

            loc_array[-1] = lend
            int_list[-1] = bend
            dtrace = algc.mag(np.concatenate([np.diff(loc_array, axis=0), np.array([[0, 0, 0]])]))

            k = np.sqrt(int_list[0] - int_list) * dtrace
            k = np.nansum(k)
            klist.append(k)
            blist.append(cb)

        klist = np.array(klist)
        blist = np.array(blist)

        return blist, klist

    def __all_b_mirror__(self):
        return True

    def __max_k__(self):
        assert (self.__all_k__()), "K not calculatable"
        return np.max(self.K)

    def __max_b__(self):
        assert (self.__all_b_mirror__()), "B not available"
        return np.nanmax(self.B)

    def __k_dir_loc__(self, loc_data, b_data, idx):

        bm = self.m_trace_b_mirror[idx]
        # print ("B: {}".format(self.m_trace_b_mirror))
        # print ("IDX in kdirloc: {}".format(idx))

        try:
            bidx = np.argmax(b_data[idx+1:] > bm) + idx + 1
        except ValueError:
            return np.NaN, bm

        int_list = b_data[idx:(bidx+1)]
        loc_list = loc_data[idx:(bidx+1)]

        den = int_list[-1] - int_list[-2]

        if den == 0:
            lw = 0
        else:
            lw = (bm - int_list[-2])/den

        bend = int_list[-2] + ((int_list[-1] - int_list[-2]) * lw)
        lend = loc_list[-2] + ((loc_list[-1] - loc_list[-2]) * lw)

        loc_list[-1] = lend
        int_list[-1] = bend

        ds = algc.mag(np.concatenate([np.diff(loc_list, axis=0), np.array([[0,0,0]])]))
        k = np.sqrt(int_list[0] - int_list) * ds
        k = np.nansum(k)

        return k, bm

    def __get_kb_for_start__(self, alpha=np.pi/2):
        if self.m_trace is not None and len(self.m_trace) > 1:

            b_m = self.B[self.start_idx]/np.sin(alpha)**2
            k = self.get_k(b_m)

            return k, b_m
        else:
            return np.nan, np.nan


    def __get_k_b_for_fl_trace_pt__(self, idx):

        loc_data, b_data = self.m_trace.copy(), self.m_trace_b_mirror.copy()
        kf, bf = self.__k_dir_loc__(idx=idx, loc_data=loc_data, b_data=b_data)

        bm = b_data[idx]

        loc_data = np.flipud(loc_data)
        b_data = np.flipud(b_data)
        bidx = len(b_data) - idx - 1

        bm2 = b_data[bidx]

        assert bm == bm2, "Wrong bidx"

        kb, bb = self.__k_dir_loc__(idx=bidx, loc_data=loc_data, b_data=b_data)

        # print ("KB: {}".format(kb))
        # print ("BB: {}".format(bb))

        k = np.nanmax([kb,kf])

        # print ("K: {}: {}".format(idx, k))
        # print ("B: {}: {}".format(idx, bf))


        return k, bf


