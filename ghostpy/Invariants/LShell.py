# File: LShell.py
# Author: Joshua Murphy
#

import numpy as np

import ghostpy.Invariants.FieldLine as fl
import ghostpy.algorithms.common as algc
import ghostpy.data.DipoleData as dpd
import ghostpy.data.GpData as gpd
import ghostpy.grids.polarGrid as pg
import scipy.interpolate as interp
import scipy.interpolate as sint
# import algorithms.FieldTracers as fts
# import algorithms.convert as algcv
# import scipy.interpolate as sint


# FIXME: Having an infinate loop is not helpfull when trying to find the line.. Need to handle invalid lines in the convergance method.
class LShell:
    def __init__(self, k=None, b_mirror=None, start_loc=None, alpha=None, data=None, save_lines=False, error_tol=1e-6, pre_converge=True):
        """
        Creates an LShell object instance, allowing several different LStar manipulations.
        :param k: Starting K for the calculation.  Only required if needing a k/loc starting point or a k/b search
        :param b_mirror: Starting B_mirror for the particle.  Only required for a B/loc starting point for a k/b search
        :param start_loc: staring location for the initial trace.  Required unless performing a k/b search
        :param data: Required for all modes. Must be of type data.GpData or one of its subclasses
        :param save_lines: True=Retain all lines for plotting. False=Discard lines, only keep trajectory
        :param error_tol: changes the error tollerance for several algorithms. Default is 1e-6
        :param pre_converge: True = calculate close approximation with 4 field line traces immediately
                             False = Only calculate the first field line.  L* may be invalid for any but Dipole Fields
                                     Until the LShell.converge_lstar() is called.
                                     This allows for more fine grained manipulation, allowing for adding of individual
                                     lines for experimentation.
        """
        np.seterr(invalid='ignore')

        assert ((k is not None and start_loc is not None) or
                (b_mirror is not None and start_loc is not None) or
                (b_mirror is None and k is None and start_loc is not None)), \
            "LShell was not called properly.\n\n " \
            "LShell must be called with one of the following three forms:\n\n" \
            "LShell(start_loc=[x,y,z], k=?, data=gp.data)\n" \
            "LShell(start_loc=[x,y,z], b_mirror=?, data=dp.data)\n" \
            "LShell(start_loc=[x,y,z], data=dp.data)\n\n" \
            "Specifiying only a starting location will result in utilizing the location B value for b_mirror.\n\n"

        assert isinstance(data, gpd.data), "The supplied data structure is of an incorrect type.  " \
                                           "All data structures must inherit ghostpy.data.GpData"


        self.k = k
        self.b = b_mirror
        self.start_l = None
        self.data = data
        self.BL = self.data.get_bisection_model()
        self.calc_boundary = self.data.get_calc_boundary()
        self.start = start_loc
        self.error_tol = error_tol
        self.valid = True
        # path tracking
        self.path = dict()
        self.conv_path = []
        self.dipole = dpd.DipoleData()
        self.pre_converge=pre_converge

        if alpha is None:
            self.alpha = np.pi/2
        else:
            self.alpha = alpha

        # line tracking
        self.save_lines = save_lines
        if self.save_lines:
            self.lines = dict()
        else:
            self.lines = None

        if self.b is not None and start_loc is not None:
            # Initialize with k/b trace at phi = 0
            assert False, "Not yet implemented"

        elif start_loc is not None and self.b is None and self.k is None:
            print("\n\n\nVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV")
            print("Calculating L*, K for point {}".format(self.start))

            newFL = self.__trace_from_location__(loc=start_loc)
            # print ("Trace Data: {}".format(newFL.trace_data_n))

            if newFL.valid and len(newFL.trace_data_n) > 0:

                self.k, self.b = newFL.__get_kb_for_start__(alpha=self.alpha)

                if self.k is None:
                    self.valid = False
                    return

                print("B_mirror = {}".format(self.b))
                print("K = {}".format(self.k))
            else:
                print ("Invalid line initiation")
                print ("Starting Location: {}".format(newFL.start))
                self.valid = False
                return

        else:
            # Initialize with location/k and find K
            newFL = self.__trace_from_location__(loc=start_loc)
            assert isinstance(newFL, fl.FieldLine), "ERROR: Incorrect type for field line."
            self.b = newFL.get_b(k=self.k)
            if self.b < 0:
                # newFL.valid = False
                self.valid = False

        if self.b > 0 and self.k >= 0:

            fp = newFL.get_footprint_rlp(re=self.calc_boundary)
            if fp is not None:
                self.path[algc.cart_to_sphere(self.start)[2]] = fp

        else:
            print("vvvvvvvvvvvvvvvvvvvvvvvv")
            print("INITIALIZATION OF LINE")
            print("Origin: {}".format(self.start))
            print("B: {}".format(self.b))
            print("K: {}".format(self.k))

            if self.k < 0:
                newFL.__all_k__(debug=True)
                print("INVALID K")
                print("NewLine all K: {}".format(newFL.K))
                print("B_Mirror: {}".format(newFL.__get_b_mirror__(self.k)))
                print("B_North: {}".format(algc.mag(newFL.trace_data_n[0:10])))
                print("B_South: {}".format(algc.mag(newFL.trace_data_s[0:10])))
                print("RE: {}".format(algc.mag(newFL.trace_n[0:10])))
            self.valid = False
            print("^^^^^^^^^^^^^^^^^^^^^^^^")

        # if we are saving lines, save it.
        if self.lines is not None and self.valid:
            self.lines[algc.cart_to_sphere(self.start)[2]] = [newFL, newFL]

        start_phi = algc.cart_to_sphere(self.start)[2]
        print ("Start Phi: {}".format(start_phi))
        self.conv_path.append(start_phi)

        # initial convergence
        if self.pre_converge:
            self.valid = self.__4_line_conv__()

        print ("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

    def __str__(self):
        return "GHOSTpy LShell object"

    def __repr__(self):
        return "-------------------\n" \
               "LShell Object:\n" \
               "Data Trace Boundary: {}\n" \
               "Data Calculation Boundary: {}\n" \
               "Starting Location: {}\n"\
               "Start RE: {}\n" \
               "- - - - - -\n"\
               "Drift Path (r, lambda, phi):\n{}\n" \
               "Number of Lines: {}\n" \
               "K: {}\n" \
               "B_mirror: {}\n" \
               "L*: {}\n" \
               "Is Valid?: {}\n" \
               "Retaining Lines: {}\n" \
               "-------------------".format(self.data.get_trace_boundary(),
                                            self.data.get_calc_boundary(),
                                            self.start,
                                            algc.mag(self.start),
                                            self.__ordered_path__(),
                                            len(self.path),
                                            self.k,
                                            self.b,
                                            self.l_star(res=1000),
                                            self.valid,
                                            isinstance(self.lines, dict))

    def l_star(self, res=1000):
        if self.valid or (not self.valid and len(self.path) >= 4):
            return 2*np.pi/self.__sub_flux__(res=res)
        else:
            return np.NaN

    def get_conv_path(self):
        return np.array(sorted(self.conv_path))

    def __sub_flux__(self, res=1000):
        """
        Utilizes Stoke's Theorem to solve the flux integral
        :param res: The number of divisions to use in the integration
        :return:
        """
        traj = np.array(self.get_drift_trajectory(res=res))
        r, lat, phi = np.hsplit(traj,3)

        r = r.flatten()
        lat = lat.flatten()
        phi = phi.flatten()

        phin = phi[0]
        phin += (np.pi*2)
        phi = np.concatenate([phi, [phin]])
        dphi = np.diff(phi)
        cosl = np.cos(lat)**2
        flux1 = (cosl * dphi)
        rflux1 = flux1/r
        flux = rflux1.sum()
        return flux

    def get_l_star_from_polar_cap_integration(self):
        return self.L_star_from_flux(flux=self.__flux_pc__(), RE=self.calc_boundary)

    def __flux_pc__(self, res=[100,1000]):
        pg1 = pg.polarGrid(RE=self.calc_boundary, resolution=res, base_lat_fun=self.__flux_pc__fun__)
        return pg1.get_surface_integral(value_fun=self.__flux_pc_value_fun)

    def __flux_pc_value_fun(self, xr, yl, zp):
        dpX, dpY, dpZ = self.dipole.get_xyz([xr, yl, zp])

        b = np.vstack(([dpX.T], [dpY.T], [dpZ.T])).T
        n = np.vstack(([xr.T], [yl.T], [zp.T])).T

        # print("B Values: {}".format(b))
        # print("N Values: {}".format(n))

        val = np.zeros_like(dpX)
        shape = np.shape(dpX)
        for x in range(shape[0]):
            for y in range(shape[1]):
                for z in range(shape[2]):
                    val[x, y, z] = b[x, y, z].dot(n[x, y, z])

        return val

    def __flux_pc__fun__(self, longitude):
        return self.__lambda_at_phi__(phi=longitude)

    def add_trace(self, phi=None):
        new_line_u, new_line_l = self.__find_line_from_phi__(phi=phi)

        if not isinstance(new_line_l, fl.FieldLine) or not isinstance(new_line_u, fl.FieldLine):
            print ("Low: {}".format(new_line_l))
            print ("High: {}".format(new_line_u))
            return False, -999

        # Add Trace Failure Modes
        if not new_line_u.valid or not new_line_l.valid:
            # TODO: work on sorting out errors from this failure
            print ("Bad Trace.  Number of Lines: {}".format(self.get_number_of_traces()))
            print ("Reason for Invalid: {}".format(new_line_u.valid_code))
            if new_line_u.valid_code == -1:
                print ("min B: {}".format(np.min(new_line_u.B)))
                print ("max B: {}".format(np.max(new_line_u.B)))
                print ("B/K response: {}".format(new_line_u.get_b(k=self.k)))
                print ("Needed B: {}".format(self.b))
                assert False, -1
            return False, new_line_u.valid_code

        # let us check the gaps, and interpolate a footprint if necessary
        if new_line_l is not None and new_line_u is not None:
            assert isinstance(new_line_u, fl.FieldLine)
            gapl = new_line_l.__get_min_b_gap__(k=self.k, b=self.b)
            gapu = new_line_u.__get_min_b_gap__(k=self.k, b=self.b)

            fpl = new_line_l.get_footprint_rlp(re=self.calc_boundary)
            fpu = new_line_u.get_footprint_rlp(re=self.calc_boundary)

            # print ("gapl: {}".format(gapl))
            # print ("gapu: {}".format(gapu))

            if gapl is not None and gapu is not None:
                if gapl - gapu == 0:
                    lw = 0
                else:
                    lw = (0 - gapu) / (gapl - gapu)
            else:
                return False, -3

            fp = fpu + (fpl-fpu)*lw

        else:
            return False, -4

        if np.any(np.isnan(fp)):
            # TODO: See if this failure still occurs
            print ("Bad Foot Print")
            return False, -5

        # If we didn't fail, we got here! Valid Line!
        self.path[phi] = fp
        print ("...")

        if self.lines is not None and new_line_u.valid:
            self.lines[phi] = [new_line_l, new_line_u]

        self.conv_path.append(phi)
        return True, 0

    def get_number_of_traces(self):
        return len(self.path)

    def get_drift_trajectory(self, res=100):
        assert res > 1, "Trajectory Resolution must be greater than 1."
        phis = np.linspace(start=0, stop=2*np.pi, num=res-1, endpoint=False)
        traj = []
        for phi in phis:
            rlam = self.__lambda_at_phi__(phi=phi, get_r=True)

            # print ("rlam {}".format(rlam))

            lam = rlam[1]
            r = rlam[0]
            # print ("Returned {} for phi = {}".format(lam, phi))
            assert lam is not None, ("phi = {} returned None.".format(phi))
            traj.append(np.array([r, lam, phi]))

        return np.array(traj)

    def get_drift_trajectory_cart(self, res=100):
        assert res > 1, "Trajectory Resolution must be greater than 1."
        if len(self.lines) == 0:
            return None
        phis = np.linspace(start=0, stop=2*np.pi, num=res-1, endpoint=False)
        traj = []
        for phi in phis:
            rlam = self.__lambda_at_phi__(phi=phi, get_r=True)
            # print ("rlam: {} for phi: {}".format(rlam, phi))
            lam = rlam[1]
            r = rlam[0]
            # print ("Returned {} for phi = {}".format(lam, phi))
            assert lam is not None, ("phi = {} returned None.".format(phi))
            traj.append(np.array(algc.sphere_to_cart(r=r, lam=lam, phi=phi)))

        return np.array(traj)

    def get_raw_path(self):
        return self.__ordered_path__()

    def converge_lstar(self, tol=1e-4, int_res=1000):
        if not self.valid:
            # TODO: Make sure we are actually invalidating only invalid lines
            print("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
            print("Bad Valid Flag in Convergence.")
            print("Number of Lines: {}".format(self.get_number_of_traces()))
            print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
            return False
        olstar_list=[]
        old_lstar = 0
        lines = np.array([0, np.pi / 2, np.pi, np.pi * 3 / 2])
        # lines = np.linspace(start=0.0, stop=2*np.pi, num=10, endpoint=False)

        new_lstar = self.l_star(res=int_res)
        olstar_list.append(old_lstar)

        count = 1
        interval = 2
        while not np.isclose(old_lstar, new_lstar, rtol=0, atol=tol, equal_nan=True):
        # while not np.any(np.isclose(olstar_list, new_lstar, rtol=0, atol=tol, equal_nan=True)):
            new_points = []

            # olstar_list.append(new_lstar)
            offset = self.__get_new_conv_offset__(pos=count, maxdiv=interval)

            new_points.append(lines[0] + offset)
            new_points.append(lines[0] + (np.pi * 2) - offset)
            new_points.append(new_points[0] + np.pi)
            new_points.append(new_points[1] - np.pi)

            for phi in new_points:
                # print ("Adding Trace: phi = {} rad".format(phi))
                self.add_trace(phi=phi)

            old_lstar = new_lstar
            new_lstar = self.l_star(res=int_res)
            print ("Number of Lines: {}".format(self.get_number_of_traces()))
            print ("New L*: {}".format(new_lstar))
            count += 1
            if count > interval:
                interval += 1
                count = 1

        return self.valid

    def converge_p2(self, depth=3):
        new_lstar = self.l_star(res=1000)

        if not self.valid:
            return False

        for x in range(int(depth)):
            phi_list = self.get_conv_list()
            for phi in phi_list:
                status, code = self.add_trace(phi=phi)
                if not status:
                    print ("Line Not Added: Code: {}".format(code))
            new_lstar = self.l_star(res=1000)
            print ("New Lstar: {}".format(new_lstar))

        return True

    def get_conv_list(self):
        phi_pre = self.get_conv_path()
        if not np.any(phi_pre ==0.0):
            phi_pre = np.concatenate([[0.0], phi_pre])
        phi_pre = np.concatenate([phi_pre, [np.pi*2]])
        new_phi = phi_pre[:-1] + (np.diff(phi_pre)/2)
        return new_phi

    @staticmethod
    def __get_new_conv_offset__(pos=0, maxdiv=10):

        assert pos != 0, "pos cannot be 0"
        if pos == maxdiv:
            maxdiv += 0.5
        den = float(pos) / float(maxdiv)
        offset = (np.pi / 2) * den

        return offset

    def __4_line_conv__(self):
        print ("Dipole L: {}".format(self.l_star(res=100)))
        line_fail = False
        lines = np.array([0, np.pi/2, np.pi, np.pi*3/2])
        for phi in lines:
            # print ("Adding Line at phi: {}".format(phi))
            # print ("Curent Phi Path: {}".format(self.path.keys()))
            line_fail = False

            if np.any(self.path.keys() == phi):
                continue
            status, code = self.add_trace(phi=phi)
            if not status:
                print("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
                print ("Line Failed to Add in 4 line convergence")
                print ("Failure Mode: {}".format(code))
                print ("Phi: {}".format(phi))
                print ("Est. L: {}".format(self.l_star(res=100)))
                print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
                line_fail = True
                self.valid = False
                break

        if not line_fail:
            return True
        else:
            return False

    @staticmethod
    def normalize_phi(phi):
        assert phi is not None, "Phi cannot be none"
        phi_1 = phi
        if phi_1 > np.pi:
            phi_1 = phi_1 - 2*np.pi
        return phi_1

    def __search_B_adaptive__(self,outerRE=30, innerRE=None, phi_0=0, ref_r=None, debug=False, tol=1e-4):

        ol = algc.mag(self.start)
        # ol = self.l(res=10)
        ib = self.data.get_trace_boundary()

        dob = 1.85 * ol
        dib = 0.65 * ol

        if dib < ib:
            dib = ib


        print ("================")
        # initial values
        phi_in = phi_out = phi_c = self.normalize_phi(phi_0)
        eps = np.finfo(np.float32).eps

        newfp = self.get_raw_path()


        lam_0 = 0
        #
        # if len(newfp) > 0:
        #     lam_0 = 0.75 * newfp[0][1]
        # else:
        #     return None, None

        print (lam_0)

        if ref_r is None:
            r_0 = r_c = self.l_star(res=100)
        else:
            r_0 = ref_r

        if debug is not None:
            stats = {}

            stats['tau_in'] = []
            stats['tau_out'] = []
            stats['tau'] = []
            stats['loc_t'] = []

        r_in = dib
        r_out = dob

        tau_in = None
        tau_out = None

        gap_out = None
        gap_in = None

        r_gap = r_out - r_in
        initcount = 0
        while tau_in is None or tau_out is None and r_gap > eps:
            # print ("Main Bounds Loop")

            loc_t = algc.sphere_to_cart(r=r_c, lam=lam_0, phi=phi_c)
            tau = self.__trace_from_location__(loc_t)

            # if we don't have a returning trace, try again.
            if len(tau.trace_n) == 0 or not algc.mag(tau.trace_n[-1]) < tau.safe_boundry:
                if tau_out is not None:
                    r_in = r_c
                else:
                    r_out = r_c

                r_c = (r_in + r_out)/2
                initcount += 1
                if initcount > 90:
                    print ("Out on init Count 1")
                    break
                continue

            # check the footprint
            tau_fp = tau.get_footprint_rlp(re=self.calc_boundary)

            # check valid footprint
            if tau_fp is None or np.any(np.isnan(tau_fp)):
                # print ("Foot Print: {}".format(tau_fp))
                # print ("Find out why a valid line is not giving a Foot print")
                return None, None
                # assert False

            # check phi_tau from both directions
            phi_gap, phi_tau = self.get_phi_gap(phi_0, tau_fp)

            # print("Phi Gap: {}".format(phi_gap))

            if False: #not np.isclose(phi_gap, 0.0, rtol=0, atol=tol):
                if phi_gap > 0:
                    phi_g_u = self.normalize_phi(phi_tau + phi_gap)
                    phi_g_l = self.normalize_phi(phi_tau - phi_gap)

                else:
                    phi_g_l = self.normalize_phi(phi_tau + phi_gap)
                    phi_g_u = self.normalize_phi(phi_tau - phi_gap)

                pcount = 0
                cthresh = 90
                while not np.isclose(phi_gap, 0.0, rtol=0, atol=tol) and pcount < cthresh:
                    # print ("Phi Loop 1")

                    phi_c = self.normalize_phi((phi_g_l + phi_g_u) / 2)

                    loc_t = algc.sphere_to_cart(r=r_c, lam=lam_0, phi=phi_c)

                    tau = self.__trace_from_location__(loc=loc_t)

                    tau_fp = tau.get_footprint_rlp(re=self.calc_boundary)
                    phi_gap, phi_tau = self.get_phi_gap(phi_0, tau_fp)

                    if phi_gap < 0:
                        phi_g_l += phi_gap
                        phi_g_l = self.normalize_phi(phi_g_l)
                    else:
                        phi_g_u += phi_gap
                        phi_g_u = self.normalize_phi(phi_g_u)

                    pcount += 1
                    if pcount > cthresh:
                        print ("Maximum Itterations on Initial Phi Search")
                        break
            else:

                b_gap_tau = tau.__get_min_b_gap__(k=self.k, b=self.b)

                # print ("B Gap Initial: {}".format(b_gap_tau))
                if b_gap_tau is not None and np.isclose(b_gap_tau, 0.0, rtol=0, atol=self.error_tol):
                    print ("Out with match.")
                    return tau, tau

                if b_gap_tau is None:
                    if gap_out is not None :
                        # print ("No gap, moving back out")
                        r_in = r_c
                    elif gap_in is None and gap_out is None:
                        # print ("Moving outer boundary in, searching for a valid point")
                        r_out = r_c
                    else:
                        r_out = r_c


                if b_gap_tau > 0:
                    if r_out >= dob:
                        dob *= 1.5
                        r_out = dob
                        r_in = 0.95 * r_c
                        if r_in < ib:
                            r_in = ib
                        r_c = (r_in + r_out)/2
                        print ("Moving outer boundary out")
                        continue
                    r_in = r_c
                    tau_in = tau
                    gap_in = b_gap_tau

                elif b_gap_tau is not None:
                    if r_in <= dib:
                        dib *= 0.5
                        r_in = dib
                        r_out = 1.1 * r_c
                        if r_in < ib:
                            r_in = ib
                        r_c = (r_in + r_out)/2
                        print ("Moving inner boundary in")
                        continue
                    r_out = r_c
                    tau_out = tau
                    gap_out = b_gap_tau

                r_c = (r_in + r_out)/2
                r_gap = r_out - r_in
                # print("R Gap: {}".format(r_gap))
            initcount += 1
            # print ("Bounds Loop Count: {}".format(initcount))
            if r_gap < eps:
                print ("Out on EPS convergence")
                return None, None
                # break

        r_out = r_out #* 2.5
        r_gap = r_out - r_in
        conv_count = 0
        nan_count = 0
        while r_gap > eps:
            # print ("Main Search Loop")
            r_c = (r_out + r_in)/2
            loc_t = algc.sphere_to_cart(r=r_c, lam=lam_0, phi=phi_c)
            tau = self.__trace_from_location__(loc=loc_t)
            tau_fp = tau.get_footprint_rlp(re=self.calc_boundary)
            phi_gap, phi_tau = self.get_phi_gap(phi_0, tau_fp)

            if False: # not np.isclose(phi_gap, 0.0, rtol=0, atol=tol):
                if phi_gap > 0:
                    phi_g_u = self.normalize_phi(phi_tau + phi_gap)
                    phi_g_l = self.normalize_phi(phi_tau - phi_gap)

                else:
                    phi_g_l = self.normalize_phi(phi_tau + phi_gap)
                    phi_g_u = self.normalize_phi(phi_tau - phi_gap)

                pcount = 0
                while not np.isclose(phi_gap, 0.0, rtol=0, atol=tol):

                    phi_c = self.normalize_phi((phi_g_l + phi_g_u) / 2)
                    loc_t = algc.sphere_to_cart(r=r_c, lam=lam_0, phi=phi_c)
                    tau = self.__trace_from_location__(loc=loc_t)
                    tau_fp = tau.get_footprint_rlp(re=self.calc_boundary)
                    phi_gap, phi_tau = self.get_phi_gap(phi_0, tau_fp)

                    if phi_gap < 0:
                        phi_g_l += (0.6*phi_gap)
                        phi_g_l = self.normalize_phi(phi_g_l)
                    else:
                        phi_g_u += (0.6*phi_gap)
                        phi_g_u = self.normalize_phi(phi_g_u)

                    pcount += 1
                    if pcount > 90:
                        # print ("Out on pcount")
                        break

            b_gap_tau = tau.__get_min_b_gap__(k=self.k, b=self.b)
            # print ("sB_gap: {}".format(b_gap_tau))

            if b_gap_tau is not None and np.isclose(b_gap_tau, 0.0, atol=self.error_tol):
                # print ("Returning on good find")
                # print("Relative B Error: {}".format(b_gap_tau))
                return tau, tau

            if b_gap_tau is None or np.isnan(b_gap_tau):
                print ("Bad Tau. Moving in.")
                r_out = (r_in + r_out)/2
                conv_count += 1
                nan_count += 1
                if nan_count > 50:
                    # print ("NaN Count: {}".format(nan_count))
                    break
                continue

            if b_gap_tau >= 0:
                # print("Looking at Tau > 0")
                # print("gap_in: {}".format(gap_in))
                # print("gap_out: {}".format(gap_out))
                if b_gap_tau < gap_in:
                    # print ("Saving Tau in with Gap: {}".format(b_gap_tau))
                    tau_in = tau
                    gap_in = b_gap_tau

                    if debug:
                        stats['tau'].append(tau)
                else:
                    pass
                    # print ("Tau: {}".format(b_gap_tau))
                r_in = r_c
            else:
                # print ("Looking at Tau < 0")
                # print("gap_in: {}".format(gap_in))
                # print("gap_out: {}".format(gap_out))

                if b_gap_tau > gap_out:
                    # print ("Saving Tau Out with Gap: {}".format(b_gap_tau))
                    tau_out = tau
                    gap_out = b_gap_tau

                    if debug:
                        stats['tau'].append(tau)
                else:
                    pass
                    # print ("Tau: {}".format(b_gap_tau))
                r_out = r_c
            r_gap = r_out - r_in
            conv_count += 1
            if debug:
                stats['tau_in'].append(tau_in)
                stats['tau_out'].append(tau_out)
            if conv_count > 150:
                print ("Out on conv count")
                break

        if debug:
            return stats
        print("Relative B Error (IN): {}".format(gap_in))
        print("Relative B Error (OUT): {}".format(gap_out))
        return tau_in, tau_out

    def get_phi_gap(self, phi_0, tau_fp):
        phi_tau = self.normalize_phi(tau_fp[2])
        phi_gap1 = self.normalize_phi(phi_0 - phi_tau)
        phi_gap2 = self.normalize_phi(phi_0 + (2 * np.pi) - phi_tau)
        phi_g_l = [phi_gap1, phi_gap2]
        phi_g_l_m = np.abs(phi_g_l)
        min_phi_loc = np.argmin(phi_g_l_m)
        phi_gap = self.normalize_phi(phi_g_l[min_phi_loc])
        return phi_gap, phi_tau

    def __search_B__(self, outerRE=15.0, innerRE=None, phi=0, direction=1):
        count = 0
        eps = np.finfo(np.float32).eps # * 10
        start_re = self.l_star(res=100)
        start_loc = algc.sphere_to_cart(r=start_re, lam=0, phi=phi)
        newLine = self.__trace_from_location__(loc=start_loc)
        assert newLine is not None, "No Line retrieved"

        ol = self.l_star(res=10)
        ib = self.data.get_trace_boundary()
        dob = 4.5 * ol
        dib = 0.25 * ol

        if dib < ib:
            dib = ib

        # print ("Epsilon: {}".format(eps))

        # print ("Start RE: {}".format(start_re))
        # print ("Footprint: {}".format(start_line.get_footprint_rlp(re=3.0)))
        # print ("Bottom South: {}".format(algc.cart_to_sphere(start_line.trace_s[-1])))
        # print ("Bottom North: {}".format(algc.cart_to_sphere(start_line.trace_n[-1])))

        # b = start_line.get_b(k=self.k)
        newRE = 0

        if innerRE is None:
            innerRE = self.data.get_trace_boundary()


        innerRE = dib
        outerRE = dob

        innerRE = np.float32(innerRE)
        outerRE = np.float32(outerRE)

        orig_outerRE = outerRE
        orig_innerRE = innerRE

        outerLine = None
        innerLine = None
        outerGap = None
        innerGap = None

        post_found = False

        # print ("InnerRE: {}".format(innerRE))
        # print ("OuterRE: {}".format(outerRE))
        search_inner = False
        gap = None
        lw_old = 0
        lw = -1

        while outerRE - innerRE > eps:
            lw_old = lw
            # print ("LW OLD: {}".format(lw_old))
            count += 1
            # print ("Loop: {}".format(count))
            if search_inner:
                newRE += 0.05
            else:

                if outerGap is not None and innerGap is not None:
                    ob = outerGap + self.b
                    ib = innerGap + self.b

                    # try:
                    #     ow = self.BL(ob)
                    #     iw = self.BL(ib)
                    #     cw= self.BL(self.b)
                    #
                    # except ValueError:
                    ow = 1.0
                    iw = 0.0
                    cw = 0.5

                    # print ("OW: {}".format(ow))
                    # print ("CW: {}".format(cw))
                    # print ("IW: {}".format(iw))
                    safe = 1.0
                    lw = ((cw-iw)/(ow-iw))*safe

                    if np.isclose(lw, lw_old, atol=5e-2,rtol=0):
                        lw = 0.5

                    if lw > 1:
                        lw = 0.9
                    #
                    # if count > 28:
                    #     lw = 0.5
                    # elif lw < 0.1:
                    #     lw = 0.1
                    # print ("LW: {}".format(lw))
                    newRE = innerRE + (outerRE-innerRE)*lw
                else:
                    newRE = (outerRE + innerRE) / 2

            # print ("newRE: {}".format(newRE))
            # print("GAP: {}".format(gap))

            if newRE < self.data.get_trace_boundary() and not search_inner:
                newRE = self.data.get_trace_boundary()*1.05
                search_inner = True

            newLoc = algc.sphere_to_cart(r=newRE, lam=0, phi=phi)
            newLine = self.__trace_from_location__(loc=newLoc)
            # newB = newLine.get_b(k=self.k)
            # print ("New B: {}".format(newB))
            #
            # if newLine.m_trace is not None:
            #     trace_m_idx = np.nanargmax(newLine.m_trace_re)
            #     trace_coord = newLine.m_trace[trace_m_idx]
            #     print("Trace Coordinate: {}".format(trace_coord))
            #     line_trace_coord = newLine.start
            #     print ("Start Line: {}".format(line_trace_coord))

            gap = newLine.__get_min_b_gap__(k=self.k, b=self.b)

            if gap is not None and np.isclose(gap, 0.0, atol=self.error_tol, rtol=0.0):
                print ("Iterations to converge: {}".format(count))
                return newLine, newLine

            # print ("Searching...")
            if gap is None:
                line_status = newLine.valid_code
                # print ("Line Status: {}".format(line_status))
                if line_status == -2:
                    outerRE = newRE
                    newLine = None
                    continue
                else:
                    if gap < 0:
                        outerRE = newRE
                        newLine = None
                    else:
                        innerRE = newRE
                        newLine = None
                    continue

            elif gap > 0:
                innerRE = newRE
                if innerGap is None or gap < innerGap:
                    # Move out
                    innerLine = newLine
                    newLine = None
                    innerGap = gap
            else:
                outerRE = newRE
                if gap > outerGap:
                    outerLine = newLine
                    newLine = None
                    outerGap = gap

        # Reverse Search Failed to get Value
        if gap is None and np.isclose(newRE, orig_outerRE, atol=eps, rtol=0):
            print ("Line Failure: Outer Boundary for phi = {}".format(phi))
            print ("Inner Gap: {}".format(innerGap))
            print ("Outer Gap: {}".format(outerGap))
            return None, None

        # Forward Search Failed to get Value
        if gap is None and np.isclose(newRE, orig_innerRE, atol=eps, rtol=0):
            print ("Line Failure: Inner Boundary for phi = {}".format(phi))
            print ("Inner Gap: {}".format(innerGap))
            print ("Outer Gap: {}".format(outerGap))
            return None, None

        # Check to see if point was bounded
        if outerGap is None or innerGap is None:
            print ("Line Failure: No Bounds for phi = {}".format(phi))
            print ("Inner Gap: {}".format(innerGap))
            print ("Outer Gap: {}".format(outerGap))
            return None, None

        print ("Inner Gap: {}".format(innerLine.__get_min_b_gap__(k=self.k, b=self.b)))
        print ("Outer Gap: {}".format(outerLine.__get_min_b_gap__(k=self.k, b=self.b)))
        print ("Iterations to converge: {}".format(count))
        return innerLine, outerLine

    def __find_line_from_phi__(self, phi=None):

        innerline = None
        outerline = None

        innerline, outerline = self.__search_B__(phi=phi)

        # if innerline is None or outerline is None:
        # innerline, outerline = self.__search_B_adaptive__(phi_0=phi)

        return outerline, innerline



    def __trace_from_location__(self, loc=None):
        if loc is None:
            # print("First Trace")
            fl1 = fl.FieldLine(data=self.data, start=self.start, error_tol=self.error_tol)
            # print ("Field Trace Complete")
            return fl1
        else:
            return fl.FieldLine(data=self.data, start=loc, error_tol=self.error_tol)

    def __lambda_at_phi__(self, phi=None, get_r=False):
        if self.get_number_of_traces() == 0:
            if not get_r:
                return np.pi/4
            else:
                return np.pi/4, self.data.get_calc_boundary()
        if self.get_number_of_traces() == 1:
            values = self.get_raw_path().flatten()
            # print ("Values: {}".format(values))
            if not get_r:
                return values[1]
            else:
                return np.array([values[0], values[1]])
        else:
            path = self.__ordered_path__()

            f = np.array(path)[:, 1]
            x = np.array(path)[:, 2]
            r = np.array(path)[:, 0]
            v = phi
            # print ("X components: {}".format(x))
            # print ("New X: {}".format(v))
            # print ("Y Components: {}".format(f))

            if v > x[-1]:
                x[0] += (2*np.pi)
                f = np.roll(f, -1)
                x = np.roll(x, -1)
                r = np.roll(r, -1)
            elif v < x[0]:
                x[-1] -= (2*np.pi)
                f = np.roll(f, 1)
                x = np.roll(x, 1)
                r = np.roll(r, 1)

            new_lambda = algc.lin_interp(v=v, X=x, F=f)
            new_r = algc.lin_interp(v=v, X=x, F=r)
            retVal = np.array([new_r, new_lambda]).flatten()

            if not get_r:
                # print ("New Lambda: {}, get_r: {}".format(new_lambda, get_r))
                return new_lambda
            else:
                # print ("RetVal: {}".format(retVal))
                return retVal


    def __ordered_path__(self):
        keys = sorted(self.path.keys())
        path = []
        for key in keys:
            path.append(self.path[key])

        return np.array(path)

    @staticmethod
    def L_star_from_flux(self, flux, Be=-3.15e4, RE=1.0):
        Rearth = RE
        L = 2 * np.pi * Be * Rearth / flux
        return L

