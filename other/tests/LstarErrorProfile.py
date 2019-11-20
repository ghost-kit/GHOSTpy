import collections as col

import fieldLine as fl
import matplotlib.pyplot as plt
import numpy as np
import paraview.simple as pv
from mpi4py import MPI
from prototype import driftShell as ds

from ghostpy.prototype import inv_common as ih


def test_error_tol(count=3):
    divs = 600
    x0 = [6.6, 0.0, 0.0]
    k = 0.0
    lt = [0, 6, 12, 18]
    tols = np.logspace(start=-7.0, stop=0.0, num=count)
    tols = tols[::-1]
    data = {}
    for tol in tols:
        print "Processing Tolerance: {}".format(tol)
        ds1 = ds.driftShell(local_times=lt, start_line=x0, K=k, mode='analytic_K', error_tol=tol)
        l = ds1.calculate_L_star(RE=1.0, lat_divs=divs,lon_divs=divs)
        data[tol] = l

    X,Y = ih.dict_to_x_y(data)
    Y2 = (np.array(Y)-6.6)/6.6 * 100

    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
    ax1.set_title("Error Tolerance vs. L* and % Error")

    ax1.set_ylabel("$L^*$")
    ax1.set_xlabel("Error Tolerance")
    line1, = ax1.semilogx(X,Y, 'ko-', label="$L^*$")

    ax2.set_ylabel("% Error")
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel("% Error")
    line2, = ax2.semilogx(X,Y2, 'rx-', label="% Error")

    plt.legend(handles=[line1, line2], loc=1)
    plt.tight_layout()
    plt.savefig("out/error_analysis/error_tol_6.6RE_4_lines.pdf")


def test_RE_based(count=3):
    dist = np.linspace(start=2, stop=8, num=count)
    divs = 600
    k = 0.0
    lt = [0, 6, 12, 18]
    tol = 1e-5

    data = {}
    for x0 in dist:
        ds1 = ds.driftShell(local_times=lt, start_line=[x0,0.0,0.0], K=k, mode='analytic_K', error_tol=tol)
        data[x0] = ds1.calculate_L_star(lat_divs=divs, lon_divs=divs)

    X,Y = ih.dict_to_x_y(data)
    Y2 = (np.array(Y)-np.array(X))/np.array(X) * 100

    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
    ax1.set_title("Starting Distance vs. % Error")

    ax1.set_ylabel("$L^*$")
    ax1.set_xlabel("Starting Distance ($R_E$)")
    line1, = ax1.plot(X,Y, 'ko-', label="$L^*$")

    ax2.set_ylabel("% Error")
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel("% Error")
    line2, = ax2.plot(X,Y2, 'rx-', label="% Error")

    plt.legend(handles=[line1, line2], loc=2)
    plt.tight_layout()
    plt.savefig("out/error_analysis/error_RE_base_4_lines.pdf")


def test_flux_grid(count=3):
    divlist = np.linspace(start=10, stop=800, num=count, dtype=np.int)
    k = 0.0
    lt = [0, 6, 12, 18]
    tol = 1e-5
    x0 = [6.6, 0.0, 0.0]
    ds1 = ds.driftShell(local_times=lt, start_line=x0, K=k, mode='analytic_K', error_tol=tol)
    data = {}
    for divs in divlist:
        l = ds1.calculate_L_star(RE=1.0, lat_divs=divs,lon_divs=divs)
        data[divs] = l
        print "L*({}): {}".format(divs, l)

    X,Y = ih.dict_to_x_y(data)
    Y2 = np.abs((np.array(Y)-6.6))/6.6 * 100

    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
    ax1.set_title("Flux Integration Grid vs. % Error\nfor 6.6 $R_E$")

    ax1.set_ylabel("$L^*$")
    ax1.set_xlabel("Flux Grid Size (n x n)")
    line1, = ax1.plot(X,Y, 'ko-', label="$L^*$")

    ax2.set_ylabel("% Error")
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel("% Error")
    line2, = ax2.semilogy(X,Y2, 'rx-', label="% Error")

    plt.legend(handles=[line1, line2], loc=4)
    plt.tight_layout()
    plt.savefig("out/error_analysis/error_Flux_grid_base_4_lines.pdf")


def test_lt_boundary(count=3):
    tol = 1e-5
    x0 = [6.6, 0.0, 0.0]
    k = 0.0
    divs = 600

    data = {}
    local_times = []
    for rounds in range(count):
        lts = np.linspace(start=0, stop=24, num=(4+rounds), endpoint=False)
        local_times.append(lts)

    for lt in local_times:
        ds1 = ds.driftShell(local_times=lt, start_line=x0, K=k, mode='analytic_K', error_tol=tol)
        data[len(lt)] = ds1.calculate_L_star(lat_divs=divs, lon_divs=divs)

    X, Y = ih.dict_to_x_y(data)
    Y2 = (np.array(Y) - 6.6) / 6.6 * 100

    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
    ax1.set_title("Number of Field Lines vs. % Error")

    ax1.set_ylabel("$L^*$")
    ax1.set_xlabel("Number of Field Lines to Define Polar Cap")
    line1, = ax1.plot(X, Y, 'ko-', label="$L^*$")

    ax2.set_ylabel("% Error")
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel("% Error")
    line2, = ax2.plot(X, Y2, 'rx-', label="% Error")

    plt.legend(handles=[line1, line2], loc=3)
    plt.tight_layout()
    plt.savefig("out/error_analysis/error_local_time_base_4_lines.pdf")


def test_latitude():
    tol = 1e-5
    lt = np.linspace(start=0, stop=24, num=2, endpoint=False)
    ds1 = ds.driftShell(local_times=lt, start_line=[6.6, 0.0, 0.0], K=0, mode='analytic_K', error_tol=tol)

    phi = ds1.get_phi_locations()
    RE=[]
    for key in phi:
        RE = ih.mag(phi[key])
        lat = ih.rad_to_deg(ih.get_lat(phi[key]))
        print "RE: {}".format(RE)
        print "lat: {}".format(lat)


def test_lstar_dipole_analytic():
    tol = 1e-5
    divs = 1200
    k = 0
    num_tests = 10
    lt = np.linspace(start=0, stop=24, num=4, endpoint=False)
    Lstars = np.linspace(start=2, stop=10, num=num_tests, endpoint=True)
    data1 = col.OrderedDict()
    error1 = col.OrderedDict()
    relerror1 = col.OrderedDict()
    intersectL = col.OrderedDict()
    intersectCalcL = col.OrderedDict()
    intError = col.OrderedDict()
    intRelError = col.OrderedDict()

    for L in Lstars:
        ds1 = ds.driftShell(local_times=lt, start_line=[L, 0.0, 0.0], K=k, mode='analytic_K', error_tol=tol)
        data1[L] = ds1.calculate_L_star(RE=1.0, lat_divs=divs, lon_divs=divs)
        error1[L] = np.abs(L - data1[L])
        relerror1[L] = np.abs(L - data1[L])/L
        intersectL[L] = ih.analytic_dipole_line_lat_crossing(L=L, r=1)
        intersectCalcL[L] = ih.get_lat(ds1.get_new_phi_for_localtime(11.34))
        intError[L] = np.abs(intersectL[L] - intersectCalcL[L])
        intRelError[L] = np.abs(intersectL[L] - intersectCalcL[L])/intersectL[L]

    # get data
    X1, Y1 = ih.dict_to_x_y(data1)
    X2, Y2 = ih.dict_to_x_y(error1)
    X3, Y3 = ih.dict_to_x_y(relerror1)
    X4, Y4 = ih.dict_to_x_y(intersectL)
    X5, Y5 = ih.dict_to_x_y(intersectCalcL)
    X6, Y6 = ih.dict_to_x_y(intError)
    X7, Y7 = ih.dict_to_x_y(intRelError)

    # plot the results
    fig1 = plt.figure()
    fig2 = plt.figure()

    axF1_1 = fig1.add_subplot(111)
    axF1_2 = fig1.add_subplot(111, sharex=axF1_1, frameon=False)

    line1, = axF1_1.plot(X1, Y1, "k*-",label="Calculated $L^{*}$")
    line2, = axF1_2.semilogy(X2, Y2, "r.-",label="Absolute Error")
    line3, = axF1_2.semilogy(X3, Y3, "bo-",label="Relative Error")
    axF1_2.yaxis.tick_right()
    axF1_2.yaxis.set_label_position("right")

    axF1_1.set_xlabel("Analytic $L^{*}$")
    axF1_1.set_ylabel("Calculated $L^{*}$")
    axF1_2.set_ylabel("Error")
    axF1_1.set_title("$L^{*}$ Error -- Analytic vs. Calculated \n on a Dipole")

    axF2_1 = fig2.add_subplot(111)
    axF2_2 = fig2.add_subplot(111, sharex=axF2_1, frameon=False)

    line4, = axF2_1.plot(X4, Y4, "k*-", label="Analytic $\lambda$ intersection")
    line5, = axF2_1.plot(X5, Y5, "r.-", label="Calculated $\lambda$ intersection")
    line6, = axF2_2.semilogy(X6, Y6, "bo-", label="Absolute Error")
    line7, = axF2_2.semilogy(X7, Y7, "k-.", label="Relative Error")

    axF2_1.set_xlabel("Analytic $L^{*}$")
    axF2_1.set_ylabel("$\lambda$ intersection ($^{\circ}$)")
    axF2_2.set_ylabel("Error")
    axF2_1.set_title("$\lambda$ Intersection Error -- Analytic vs. Calculated \n on a Dipole")
    axF2_2.yaxis.tick_right()
    axF2_2.yaxis.set_label_position("right")

    axF1_1.legend(handles=[line1, line2, line3], loc=2, prop={'size': 8})
    fig1.tight_layout()
    fig1.savefig("out/error_analysis/Calculated_L_Error.pdf")

    axF2_1.legend(handles=[line4, line5, line6, line7], loc=2, prop={'size': 8})
    fig2.tight_layout()
    fig2.savefig("out/error_analysis/Lambda_intersection_error.pdf")


def test_grid_vs_anal():
    print "Conducting grid vs. analytic dipole test"
    # files_res = [3,5,7,9] #,11,13] #,15,17]
    files_res = [9]
    test_re = np.linspace(start=2, stop=10, num=16)
    min_error = col.OrderedDict()
    max_error = col.OrderedDict()
    med_error = col.OrderedDict()
    for re in test_re:
        print "Running on L*={}".format(re)
        int_anal = col.OrderedDict()
        int_grid = col.OrderedDict()
        int_err = col.OrderedDict()

        for res in files_res:
            print "Executing resolution {}".format(res)
            filename="../../out/fields/dipole_dp0_{}_grid.vts".format(res)
            data = pv.OpenDataFile(filename=filename)
            fl1 = fl.fieldLine(data=data, start=[re, 0.0, 0.0])
            int_grid[res] = ih.get_lat(fl1.get_location_for_RE(RE=1.0))
            int_anal[res] = ih.analytic_dipole_line_lat_crossing(L=re, r=1.0)
            int_err[res] = np.abs(int_anal[res] - int_grid[res])

        # X1,Y1 = ih.dict_to_x_y(int_grid)
        # X2,Y2 = ih.dict_to_x_y(int_anal)
        # X3,Y3 = ih.dict_to_x_y(int_err)
        #
        # fig1 = plt.figure()
        # ax1 = fig1.add_subplot(111)
        # ax2 = ax1.twinx()
        #
        # line1, = ax1.plot(X1,Y1, 'bo-', label="Grid Trace")
        # line2, = ax1.plot(X2,Y2, 'r*-', label="Analytic Solution")
        # line3, = ax2.semilogy(X3,Y3, 'k.-', label="Error")
        #
        # ax1.set_title("Convergence of $\lambda$ Intersection in a Gridded Dipole Field\nfor $^* =\ {}$".format(re))
        # ax1.set_xlabel("Grid Cells per $R_E$")
        # ax1.set_ylabel("$\lambda$ Intersection")
        # ax2.set_ylabel("Absolute Error")
        #
        # ax1.legend(handles=[line1, line2, line3], loc=4)
        # fig1.tight_layout()

        # fig1.savefig("out/error_analysis/Lambda_intersection_{}_grid.pdf".format(re))
        # plt.close(fig1)
        min_error[re] = np.min(int_err.values())
        max_error[re] = np.max(int_err.values())
        med_error[re] = np.median(int_err.values())

    #X1, Y1 = ih.dict_to_x_y(min_error)
    #X2, Y2 = ih.dict_to_x_y(max_error)
    X3, Y3 = ih.dict_to_x_y(med_error)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    #line1, = ax1.plot(X1, Y1, 'bo-', label="Min Error")
    #line2, = ax1.plot(X2, Y2, 'r*-', label="Max Error")
    line3, = ax1.plot(X3, Y3, 'k.-', label="Error")

    ax1.set_title("Error Rates for Field Line Tracing by $L^*$\nGrid Cells per $\\text{R_E = 9$ ")
    ax1.set_xlabel("$L^*$")
    ax1.set_ylabel("Absolute Error")

    # ax1.legend(handles=[line1, line2, line3], loc='best')
    fig1.tight_layout()

    fig1.savefig("out/error_analysis/Lambda_error_rates_by_re.pdf")


def test_grid_vs_anal_L(res=9):
    test_re = np.linspace(start=4.0, stop=10.0, num=20)
    med_error = col.OrderedDict()

    filename = "../../out/fields/dipole_dp0_{}_grid.vts".format(res)
    data = pv.OpenDataFile(filename=filename)
    fl1 = fl.fieldLine(data=data, start=[2.5, 0.0, 0.0])

    for re in test_re:
        print "Running on L*={}".format(re)
        fl1.recompute_field_line(new_start_location=[re, 0.0, 0.0])
        int_grid = ih.get_lat(fl1.get_location_for_RE(RE=2.5))
        int_anal = ih.analytic_dipole_line_lat_crossing(L=re, r=2.5)
        int_err = np.abs(int_anal - int_grid)
        med_error[re] = int_err

    X3, Y3 = ih.dict_to_x_y(med_error)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.semilogy(X3, Y3, 'k.-', label="Error")

    ax1.set_title(" $\lambda$ Intersection Error Rates for Field Line Tracing (by $L^*$)\n{} grid cell per $R_E$ Dipole Field".format(res))
    ax1.set_xlabel("$L^*$")
    ax1.set_ylabel("Absolute Error")

    # ax1.legend(handles=[line1, line2, line3], loc='best')
    fig1.tight_layout()
    fig1.savefig("out/error_analysis/Lambda_error_rates_by_re_rect.pdf")


def test_grid_vs_anal_L_lfm_grid(res=9):
    test_re = np.linspace(start=4.0, stop=10.0, num=20)
    med_error = col.OrderedDict()

    filename = "test_data/lfm_dipole_test.vts"
    data = pv.OpenDataFile(filename=filename)
    fl1 = fl.fieldLine(data=data, start=[2.5, 0.0, 0.0])

    for re in test_re:
        print "Running on L*={}".format(re)
        fl1.recompute_field_line(new_start_location=[re, 0.0, 0.0])
        int_grid = ih.get_lat(fl1.get_location_for_RE(RE=2.5))
        int_anal = ih.analytic_dipole_line_lat_crossing(L=re, r=2.5)
        int_err = np.abs(int_anal - int_grid)
        med_error[re] = int_err

    X3, Y3 = ih.dict_to_x_y(med_error)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.semilogy(X3, Y3, 'k.-', label="Error")

    ax1.set_title("$\lambda$ Intersection Error Rates for Field Line Tracing (by $L^*$)\nLFM True Double grid Dipole Field")
    ax1.set_xlabel("$L^*$")
    ax1.set_ylabel("Absolute Error")

    # ax1.legend(handles=[line1, line2, line3], loc='best')
    fig1.tight_layout()
    fig1.savefig("out/error_analysis/Lambda_error_rates_by_re_lfm.pdf")



comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 6:
    test_error_tol(count=10)
elif rank == 1:
    test_RE_based(count=10)
elif rank == 2:
    test_flux_grid(count=4)
elif rank == 3:
    test_lt_boundary(count=10)
elif rank == 4:
    test_latitude()
elif rank == 5:
    test_lstar_dipole_analytic()
elif rank == 0:
    test_grid_vs_anal()

print "Rank {} is Done.".format(rank)




