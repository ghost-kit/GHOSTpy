from __future__ import absolute_import

import numpy as np
import unittest as ut

import ghostpy.Invariants.LShell as ls
import ghostpy.algorithms.convert as algx
import ghostpy.algorithms.common as algc
import ghostpy.algorithms.FieldTracers as ft

import ghostpy.data.VtkData as vdt
import ghostpy.data.LfmVtkData as lvdt
import ghostpy.data.DipoleData as dpd
import ghostpy.data.GpData as gpd
import ghostpy.plotting.FieldLinePlot as flplt
import ghostpy.Invariants.FieldLine as fl
import csv


class TestInvariantUnits(ut.TestCase):

    # def test_set_start_location(self):
    #     print("Testing set_start_location")
    #     field_line = fl.FieldLine(start=[6.6, 0, 0], data=dpd.DipoleData())
    #     self.assertItemsEqual(field_line.get_start_location(), [6.6, 0, 0])
    #     field_line.set_start_location([6.6, 1.0, 3.3])
    #     self.assertItemsEqual(field_line.get_start_location(), [6.6, 1.0, 3.3])
    #     print("Set Start Location Test Complete")
    #
    # def test_get_trace_mode(self):
    #     print("Testing get_trace_mode")
    #     field_line = fl.FieldLine(start=[6.6, 0, 0], data=dpd.DipoleData())
    #     self.assertEqual(field_line.get_trace_mode(), 'dipole')
    #     print("Get Trace Mode Test Complete")
    #
    # def test_data(self):
    #     print("Testing Data Modules")
    #     field_line = fl.FieldLine(start=[6.6, 0, 0], data=dpd.DipoleData())
    #     self.assertEqual(field_line.get_trace_mode(), 'dipole')
    #     # self.assertEqual(self.field_line_lfm.get_trace_mode(), 'vtkData')
    #     print("Data Modules Test Complete")
    #
    # def test_get_xyz(self):
    #     print("Test get_xyz from Data")
    #     data = vdt.VtkData(filename="./test_data/lfm_dipole_test.vts", vector="B")
    #     data2 = dpd.DipoleData()
    #
    #     val = data.get_xyz((6.6, 0.0, 0.0))
    #     val2 = data2.get_xyz((6.6, 0.0, 0.0))
    #
    #     # print ("Gridded: {}\nAnalytic: {}".format(val, val2))
    #     val_m = algc.mag(val)
    #     val2_m = algc.mag(val2)
    #
    #     self.assertGreater(algc.mag(val), 0.0, msg="Dipole Data Failed")
    #     self.assertGreater(algc.mag(val2), 0.0, msg="Dipole Analytic Data Failed")
    #     self.assertAlmostEqual(algc.mag(val), algc.mag(val2), 0, msg="Values Not Almost Equal. val1: {}, val2: {}".format(val_m, val2_m))
    #     print("Test get_xyz Complete")
    #
    # def test_min_b(self):
    #     print("Testing __min_B__ routine")
    #     field_line = fl.FieldLine(start=[6.6, 0, 0], data=dpd.DipoleData())
    #     min_n, min_s = field_line.__min_B__()
    #     print("__min_B__ routine Test Complete")
    #
    # def test_resample_line(self):
    #     print("Testing Line Resample")
    #     fl2 = fl.FieldLine(start=[6.6, 1.6, 3.5], data=dpd.DipoleData(), error_tol=1e-4)
    #     fl2.__resample_line__()
    #
    #     post_n = fl2.trace_data_n
    #     post_s = fl2.trace_data_s
    #
    #     post_mag_n = algc.mag(post_n)
    #     post_mag_s = algc.mag(post_s)
    #     # print ("North Values Monotonically increasing: {} (ignoring last value - may be undefined)".format(algc.mon_inc(post_mag_n[:-1])))
    #     # print ("South Values Monotonically increasing: {} (ignoring last value - may be undefined)".format(algc.mon_inc(post_mag_s[:-1])))
    #
    #     # TODO: Fix so tolorance can be much better.
    #     self.assertTrue(np.all(np.isclose(post_mag_s, post_mag_n, atol=1, rtol=0.0, equal_nan=True)), "North vs. South Asserted values failed:\n {}".format(np.isclose(post_mag_s, post_mag_n, atol=1e-1, rtol=0.0)))
    #
    #     fl3 = fl.FieldLine(start=[6.6, 0.0, 3.5], data=dpd.DipoleData())
    #     fl3.__resample_line__()
    #
    #     post_n = fl3.trace_data_n
    #     post_s = fl3.trace_data_s
    #
    #     # print ("North Values Monotonically increasing: {} (ignoring last value - may be undefined)".format(algc.mon_inc(post_mag_n[:-1])))
    #     # print ("South Values Monotonically increasing: {} (ignoring last value - may be undefined)".format(algc.mon_inc(post_mag_s[:-1])))
    #
    #     post_mag_n = algc.mag(post_n)
    #     post_mag_s = algc.mag(post_s)
    #     # print (post_mag_n)
    #     # print (post_mag_s)
    #
    #     self.assertTrue(np.all(np.isclose(post_mag_s[:161], post_mag_n[:161], atol=1, rtol=0.0, equal_nan=True)), "asserted values failed:\n {}".format(np.isclose(post_mag_s, post_mag_n, atol=1e-1, rtol=0.0)))
    #     # print ("Passed re-sample test for non-aligned traces. Starting position: {} and {}".format(fl2.start, fl3.start))
    #
    #     fl4 = fl.FieldLine(start=[6.6, 0.0, 0.0], data=dpd.DipoleData())
    #     fl4.__resample_line__()
    #
    #     post_n = fl4.trace_data_n
    #     post_s = fl4.trace_data_s
    #
    #     # print ("North Values Monotonically increasing: {} (ignoring last value - may be undefined)".format(algc.mon_inc(post_mag_n[:-1])))
    #     # print ("South Values Monotonically increasing: {} (ignoring last value - may be undefined)".format(algc.mon_inc(post_mag_s[:-1])))
    #
    #     post_mag_n = algc.mag(post_n)
    #     post_mag_s = algc.mag(post_s)
    #     # print (post_mag_n)
    #     # print (post_mag_s)
    #
    #     self.assertTrue(np.all(np.isclose(post_mag_s[:161], post_mag_n[:161], atol=1, rtol=0.0, equal_nan=True)), "asserted values failed:\n {}".format(np.isclose(post_mag_s, post_mag_n, atol=1e-1, rtol=0.0)))
    #     # print ("Passed re-sample test for aligned traces. Starting position: {}".format(fl4.start))
    #     # print("Line Re-sample Test complete.")
    #     print("Line Resample Test Complete")
    #
    # def test_mag(self):
    #     print("Testing magnitude algorithm")
    #     test = np.array([[1, 1, 1], [2, 2, 2], [3, 3, 3]])
    #     test2 = np.array([1, 1, 1])
    #     mags = algc.mag(test)
    #     mag2 = algc.mag(test2)
    #     self.assertEqual(mag2, mags[0], "Failure of mode selection within magnitude algorithm")
    #     print("Magnitude Algorithm Test complete")

    # def test_calc_k(self):
    #     print("Testing __kb_model__ routine")
    #     data_lfm = vdt.VtkData(filename="./test_data/WHIDouble.vts", vector='B')
    #     # fl1 = fl.FieldLine(start=[8.819081929894001524,  0.0, 0.0], data=data_lfm)
    #     fl1 = fl.FieldLine(start=[0.42526826, -0.88136268,  2.73213792], data=data_lfm)
    #
    #     k = fl1.__kb_model__()
    #     K = fl1.__get_k__(106)
    #     B = fl1.__get_b_mirror__(K)
    #     print ("K: {}".format(K))
    #     print ("B: {}".format(B))
    #     self.assertTrue(k)
    #     print("__kb_model__ Test Complete")
    #
    # def test_get_b_mirror(self):
    #     print("Testing __b_mirror__ routine")
    #     fl1 = fl.FieldLine(start=[6.6, 0.0, 0.0], data=dpd.DipoleData())
    #     k0 = fl1.__b_mirror__(k=0)
    #     k330000 = fl1.__b_mirror__(k=330000)
    #     self.assertEqual(k330000, -1)
    #     self.assertGreater(k0, 0)
    #     print("B mirror test complete")
    #
    # def test_get_k(self):
    #     print("Testing __calculate_K__ Routine")
    #     fl1 = fl.FieldLine(start=[6.6, 0.0, 0.0], data=dpd.DipoleData())
    #     b100 = fl1.__get_k__(b_mirror=100)
    #     b2600 = fl1.__get_k__(b_mirror=2600)
    #     b30000 = fl1.__get_k__(b_mirror=3000)
    #     # print(fl1.__max_k__())
    #     # print(fl1.__max_b__())
    #     # print("Inner Boundary: {}".format(fl1.data.get_inner_boundary()))
    #     # print("line: {}".format(algc.mag(fl1.trace_data_n)))
    #     # print("line: {}".format(algc.mag(fl1.trace_data_s)))
    #     # print("Position: {}".format(algc.mag(fl1.trace_n)))
    #     # print("Position: {}".format(algc.mag(fl1.trace_s)))
    #     # print("K: {}".format(fl1.K))
    #     # print("B: {}".format(fl1.B))
    #     self.assertEqual(b100,-1)
    #     self.assertIsNotNone(b2600)
    #     self.assertIsNotNone(b30000)
    #     print("__calculate_K__ Routine Test Complete")
    #
    # def test_get_intersect(self):
    #     print("Testing Field line intersection routine")
    #     re = 3.0
    #     fl1 = fl.FieldLine(start=[6.6, 0.0, 0.0], data=dpd.DipoleData())
    #     intxyz = fl1.get_footprint_xyz(re=re)
    #     intrlp = fl1.get_footprint_rlp(re=re)
    #
    #     # print ("XYZ: {}\nRLP: {}".format(algc.mag(intxyz), intrlp[0]))
    #     self.assertIsNotNone(intrlp)
    #     self.assertIsNotNone(intxyz)
    #     self.assertAlmostEqual(intrlp[0], re, places=3)
    #     self.assertAlmostEqual(algc.mag(intxyz), re, places=3)
    #     self.assertGreaterEqual(intrlp[1], 0, "ERROR: We are reading the path in the wrong hemisphere!!\nLatitude Received: {} degrees".format(intrlp[1] * 180 / np.pi))
    #     print("Intersection Point (Footprint) test Complete")
    #
    # # def test_LShell_init(self):
    # #     print("Testing LShell_init")
    # #     ls1 = ls.LShell(k=0, start_loc=[6.6, 0.0, 0.0], data=dpd.DipoleData())
    # #     # print(repr(ls1))
    # #     print("LShell init test complete")
    #
    # # def test_get_number_of_traces(self):
    # #     print("Testing get_number_of_traces Method")
    # #     ls1 = ls.LShell(k=0, start_loc=[6.6, 0.0, 0.0], data=dpd.DipoleData())
    # #     num_traces = ls1.get_number_of_traces()
    # #     self.assertEquals(num_traces, 1, msg="Wrong number of traces returned. Returned {}, expected 1.".format(num_traces))
    # #     # print ("Number of Traces: {}".format(num_traces))
    # #     print("Number of Traces Test Complete")
    #
    #
    # def test_sphere_2_cart(self):
    #     phi = 11 * np.pi / 6
    #     lam = np.pi/3
    #     r = 1.0
    #     cart = algc.sphere_to_cart(r=r, lam=lam, phi=phi)
    #     sphere = algc.cart_to_sphere(cart)
    #
    #     self.assertAlmostEqual(r, sphere[0], places=12)
    #     self.assertAlmostEqual(lam, sphere[1], places=12)
    #     self.assertAlmostEqual(phi, sphere[2], places=12)
    #
    # # def test_add_trace(self):
    # #     print ("Testing add_trace Method")
    # #     # data_dpg = vdt.VtkData(filename="./test_data/lfm_t96_test.vts", vector="B")
    # #     data_dpg = dpd.DipoleData()
    # #     # ls1 = ls.LShell(k=0, start_loc=[6.6, 0.0, 0.0], data=dpd.DipoleData())
    # #     ls1 = ls.LShell(k=0, start_loc=[6.6, 0.0, 0.0], data=data_dpg)
    # #     ls1.add_trace(phi=np.pi)
    # #     ls1.add_trace(phi=np.pi/4)
    # #     ls1.add_trace(phi=np.pi/6)
    # #     ls1.add_trace(phi=np.pi*8/6)
    # #     print(repr(ls1))
    # #     print("Adding Traces Complete")
    # #
    # #     print ("add_trace testing Complete")
    #
    #
    # def test_vtk_field_trace(self):
    #     data_dpg = vdt.VtkData(filename="./test_data/lfm_t96_test.vts", vector="B")
    #     # data_dpg = vdt.VtkData(filename="./test_data/lfm_dipole_test.vts", vector="B")
    #     # data_dpg = vdt.VtkData(filename="./test_data/dipole_dp0_3_grid.vts", vector="B")
    #
    #     # data_dpg = dpd.DipoleData()
    #     ls1 = ls.LShell(k=250, start_loc=[7.1, 1.2, 0], data=data_dpg, error_tol=1e-7, save_lines=True)
    #     print(repr(ls1))
    #     keys = ls1.lines.keys()
    #     # print ("Keys: {}".format(keys))
    #     B = ls1.lines[keys[0]].B
    #     K = ls1.lines[keys[0]].K
    #
    #     # lines = np.linspace(start=0, stop=np.pi * 2, num=48, endpoint=False)
    #     lines = [np.pi / 2, np.pi, np.pi * 3 / 2, 0]
    #     for line in lines:
    #         success = ls1.add_trace(phi=line)
    #         if success:
    #             # print("L* from full polar cap integration: {}".format(ls1.get_l_star_pc()))
    #             print(repr(ls1))
    #         else:
    #             print ("Line does not exist at the provided phi: {}".format(line))
    #             print (repr(ls1))
    #
    #     # lkeys = ls1.lines.keys()
    #     flplot = flplt.FieldLinePlot()
    #     flplot.add_title("L* = {}".format(ls1.l(res=10000)))
    #
    #     flplot.plot_shell_field_lines(lshell=ls1)
    #     flplot.show()
    #
    #
    # def test_rotate_y(self):
    #     # Data = lvdt.LfmVtkData(filename="./test_data/lfm_dipole_test.vts", vector="B")
    #     x = algx.__rotate_y__([1, 0, 2], -np.pi / 2)
    #     print(x)

    #
    # def test_get_b_mirror_list(self):
    #     data_lfm = vdt.VtkData(filename="./test_data/WHIQuad.vts", vector='B')
    #     fl1 = fl.FieldLine(start=[8.7,0,0], data=data_lfm)
    #
    #     list = fl1.__get_b_mirror_list__(k=100)
    #
    #     print ("List: {}".format(list))
    #


    def test_x_axis(self):
        # data_dip = vdt.VtkData(filename="./test_data/lfm_dipole_test.vts", vector="B")
        data_lfm = vdt.VtkData(filename="./test_data/lfm_t96_test.vts", vector="B")
        # data_lfm = vdt.VtkData(filename="/Users/jjm390/src/ghostpy/ghostpy/samples/preprocess/outdata/MHD_148.vts", vector='B')

        re = np.linspace(start=-14.5, stop=14.5, num=50)
        pos = []
        for r in re:
            pos.append([r, 0, 0])

        idx = 0
        for loc in pos:
            # print("\n\n\nCalculating L* for location: {} (index {})".format(loc, idx))
            # lsd = ls.LShell(start_loc=loc, data=data_dip, save_lines=True, error_tol=1e-6)
            # lst = ls.LShell(start_loc=loc, data=data_t96, save_lines=True, error_tol=1e-6)
            lsl = ls.LShell(start_loc=loc, data=data_lfm, save_lines=True, error_tol=1e-6)

            # print ("Pre Convergence Dipole L*: {}".format(lsd.l(res=10000)))
            # print ("Pre Convergence T96 L*: {}".format(lst.l(res=10000)))
            print ("Pre Convergence LFM L*: {}".format(lsl.l_star(res=10000)))

            # flpd = flplt.FieldLinePlot()
            # flpd.add_title("Pre Convergence Dipole $L^*$ = {}, $K$= {}, $B_m$= {}\nStarting Location: {} (RE)".format(lsd.l(res=10000),  lsd.k, lsd.b, lsd.start))
            # flpd.plot_shell_field_lines(lshell=lsd)
            # flpd.plot_drift_boundary(lshell=lsd)
            # flpd.savePDF(filename="./out2/lstar_dip_sl{}_pre.pdf".format(loc))
            #
            # flpt = flplt.FieldLinePlot()
            # flpt.add_title("Pre Convergence T96 $L^*$ = {}, $K$= {}, $B_m$= {}\nStarting Location: {} (RE)".format(lst.l(res=10000),  lst.k, lst.b, lst.start))
            # flpt.plot_shell_field_lines(lshell=lst)
            # flpt.plot_drift_boundary(lshell=lst)
            # flpt.savePDF(filename="./out2/lstar_t96_sl{}_pre.pdf".format(loc))
            #
            # flpl = flplt.FieldLinePlot()
            # flpl.add_title("Pre Convergence LFM $L^*$ = {}, $K$= {}, $B_m$= {}\nStarting Location: {} (RE)".format(lsl.l(res=10000),  lsl.k, lsl.b, lsl.start))
            # flpl.plot_shell_field_lines(lshell=lsl)
            # flpl.plot_drift_boundary(lshell=lsl)
            # flpl.savePDF(filename="./out2/lstar_lfm_sl{}_pre.pdf".format(loc))

            # lsd.converge_lstar(tol=1e-4)
            # lst.converge_lstar(tol=1e-4)
            lsl.converge_lstar(tol=1e-3)
            #
            # print ("Post Convergence Dipole L*: {}".format (lsd.l(res=10000)))
            # print ("Post Convergence T96 L*: {}".format(lst.l(res=10000)))
            print ("Post Convergence LFM L*: {}".format(lsl.l_star(res=10000)))
            #
            # flpd = flplt.FieldLinePlot()
            # flpd.add_title("Post Convergence Dipole $L^*$ = {}, $K$= {}, $B_m$= {}\nStarting Location: {} (RE)".format(lsd.l(res=10000),  lsd.k, lsd.b, lsd.start))
            # flpd.plot_shell_field_lines(lshell=lsd)
            # flpd.plot_drift_boundary(lshell=lsd)
            # flpd.savePDF(filename="./out2/lstar_dip_sl{}_post.pdf".format(loc))
            #
            # flpt = flplt.FieldLinePlot()
            # flpt.add_title("Post Convergence T96 $L^*$ = {}, $K$= {}, $B_m$= {}\nStarting Location: {} (RE)".format(lst.l(res=10000),  lst.k, lst.b, lst.start))
            # flpt.plot_shell_field_lines(lshell=lst)
            # flpt.plot_drift_boundary(lshell=lst)
            # flpt.savePDF(filename="./out2/lstar_t96_sl{}_post.pdf".format(loc))
            #
            # flpl = flplt.FieldLinePlot()
            # flpl.add_title("Post Convergence LFM $L^*$ = {}, $K$= {}, $B_m$= {}\nStarting Location: {} (RE)".format(lsl.l(res=10000),  lsl.k, lsl.b, lsl.start))
            # flpl.plot_shell_field_lines(lshell=lsl)
            # flpl.plot_drift_boundary(lshell=lsl)
            # flpl.savePDF(filename="./out2/lstar_lfm_sl{}_post.pdf".format(loc))

            idx += 1


    # def test_large_area_point(self):
    #     data_lfm = vdt.VtkData(filename="./test_data/WHIDouble.vts", vector='B')
    #     points = data_lfm.points
    #
    #     RE = algc.mag(points)
    #     pointsL = RE < 3.5
    #     pointsU = RE > data_lfm.get_calc_boundary()
    #
    #     ptsu = np.logical_and(pointsL, pointsU)
    #     pts = points[np.where(ptsu)]
    #
    #     print ("Number of Points: {}".format(np.shape(pts)[0]))
    #     count = 0
    #     for loc in pts[0:10]:
    #         # print("Calculating L* for: {}".format(loc))
    #         lsl = ls.LShell(start_loc=loc, data=data_lfm, save_lines=False, error_tol=1e-6)
    #         print("[{}] L*({}): {}".format(count, loc, lsl.l(res=10000)))
    #         count += 1
    #

    # def test_conv_p2(self):
    #     data_lfm = vdt.VtkData(filename="./test_data/lfm_dipole_test_single.vts", vector='B')
    #     loc = [-1.9372226,   9.96556568,  3.6239264 ]
    #     lsl = ls.LShell(start_loc=loc, data=data_lfm)
    #
    #     print (repr(lsl))
    #
    #     lsl.converge_p2(tol=1e-5)
    #     print ("conv_path: {}".format(lsl.get_conv_path()))
    #
    #     print (repr(lsl))


    # def test_nan_point_problem(self):
    #     print ("Testing Specific Points")
    #     # data_lfm = vdt.VtkData(filename="./test_data/lfm_dipole_test_quad.vts", vector='B')
    #     # data_lfm = vdt.VtkData(filename="./test_data/WHISingle.vts", vector='B')
    #     data_lfm = vdt.VtkData(filename="./test_data/WHIDouble.vts", vector='B')
    #     # data_lfm = vdt.VtkData(filename="./test_data/lfm_dipole_test_quad.vts", vector='B')
    #
    #     pos = [[-1.60911155, 1.912817, 0.63764501] ,
    #            [1.85464418, 2.00403666, 1.17864108],
    #            [2.42400000e+00, 2.00000000e-12, -2.00000000e-12],
    #            [2.54356503, 0.49542972, 2.19501495],
    #            [-14.50, 0.0, 0.0]]
    #
    #     # pos = [[1.3328408, 0.19456816, -2.19009161]]
    #     # pos = [[7.26819780e+00, 4.19629596e+00, -2.00000000e-12]]
    #     # pos = [[-2.17350745, -0.71550012, -2.61292052]]
    #     # pos = [[-4.54661465,  0.37921321,  4.26848459]]
    #     # pos = [[-2.40793252,  0.18940264, -0.58712971]]
    #     # pos = [[-2.27353,	0,	-1.22493]]
    #     # pos = [[-0.86630422, -0.56160945,  2.48822451]]
    #     # pos = [[6.22623232e+00,   2.00000000e-12,  -2.00000000e-12]]
    #     # pos = [[-12.0,   0.0,  0]]
    #
    #
    #     for loc in pos:
    #         print("Calculating position {}".format(loc))
    #         lsl = ls.LShell(start_loc=loc, data=data_lfm, save_lines=True, error_tol=2e-6)
    #         print ("L* = {}".format(lsl.l(res=10000)))
    #         lsl.converge_lstar(tol=1e-2)
    #         # lsl.converge_p2(depth=4)
    #         print (repr(lsl))
    #
    #         flpl = flplt.FieldLinePlot()
    #         flpl.add_title("LFM $L^*$ = {}, $K$= {}, $B_m$= {}\nStarting Location: {} (RE)".format(lsl.l(res=10000),  lsl.k, lsl.b, lsl.start))
    #         flpl.plot_shell_field_lines(lshell=lsl)
    #         flpl.plot_drift_boundary(lshell=lsl)
    #         flpl.show()
    #         # flpl.savePDF(filename="./out2/lstar_lfm_sl{}.pdf".format(loc))

    # def test_convergance(self):
    #
    #     # data_dpg = vdt.VtkData(filename="./test_data/lfm_dipole_test.vts", vector="B")
    #     # data_dpg = vdt.VtkData(filename="./test_data/lfm_t96_test.vts", vector="B")
    #     data_dpg = vdt.VtkData(filename="./test_data/lfm_lfmdata_test.vts", vector='B')
    #
    #     csvfile = open("./test_data/rbsp.csv")
    #     sat_dat = csv.DictReader(csvfile)
    #     pos = []
    #     sat = []
    #     for row in sat_dat:
    #         keys = row.keys()
    #         pos.append(algx.km_to_re(np.array([float(row['X']), float(row['Y']), float(row['Z'])])))
    #         sat.append(str(row['Sat']))
    #
    #     pos = np.array(pos)
    #
    #     posb = pos[np.where(np.array(sat)=='rbspb')]
    #     posa = pos[np.where(np.array(sat)=='rbspa')]
    #
    #     posaRE = algc.mag(posa)
    #     posbRE = algc.mag(posb)
    #
    #     posaRE_out = posaRE > data_dpg.get_trace_boundary()
    #     posbRE_out = posbRE > data_dpg.get_trace_boundary()
    #
    #     posa = posa[np.where(posaRE_out)]
    #     posb = posb[np.where(posbRE_out)]
    #
    #     # print ("Position of RBSP-A: \n{}".format(posa))
    #     # print ("Position of RBSP-B: \n{}".format(posb))
    #
    #     posRE = algc.mag(pos)
    #     posOUT = posRE > data_dpg.get_trace_boundary()
    #
    #     pos = pos[np.where(posOUT)]
    #
    #     count = 0
    #     for loc in np.concatenate([posa, posb])[3:]:
    #         print("Calculating L* for: {}, index {}".format(loc, count))
    #         ls1 = ls.LShell(start_loc=loc, data=data_dpg, error_tol=1e-6, save_lines=True)
    #         print (repr(ls1))
    #         print ("Converging L")
    #         ls1.converge_lstar(tol=1e-3)
    #         print (repr(ls1))
    #         if ls1.valid:
    #             flplot = flplt.FieldLinePlot()
    #             flplot.add_title("$L^*$ = {}, $K$= {}, $B_m$= {}\nStarting Location: {} (RE)".format(ls1.l(res=10000),  ls1.k, ls1.b, ls1.start))
    #             flplot.plot_shell_field_lines(lshell=ls1)
    #             flplot.plot_drift_boundary(lshell=ls1)
    #             # flplot.ax.scatter(posb[:, 0], posb[:, 1], posb[:, 2], 'y--')
    #             # flplot.ax.scatter(posa[:, 0], posa[:, 1], posa[:, 2], 'r-.')
    #             # flplot.savePDF(filename="./out/lstar_Dipole_sl{}.pdf".format(loc))
    #             # flplot.savePDF(filename="./out/lstar_T96_sl{}.pdf".format(loc))
    #             flplot.savePDF(filename="./out/lstar_lfm_sl{}.pdf".format(loc))
    #         count +=1
    #
