# File: DriftShellErrorAnalysis.py
# Author: Joshua Murphy
# Project: PhD Dissertation
# Purpose: To identify the error rates in various parts of the drift-shell calculations (DIPOLE)
# ==============================================================================================

import collections as col
import datetime as dt

import numpy as np
import paraview.simple as pv
from mpi4py import MPI
from prototype import driftShell as ds

from ghostpy.prototype import inv_common as ih

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

lines_list = [4, 8, 16]  # 32, 64, 128, 256, 512]

for lines in lines_list:

    calcTimes = np.linspace(0, 24, num=lines, endpoint=False)

    start_location = np.linspace(6.8, 3.0, num=size, endpoint=True)
    grid_files = ["dipole_dp0_3_small_grid.vts"] #, "dipole_dp0_10_large_grid.vts"]

    plot_data = col.OrderedDict()

    for grid_file in grid_files:
        data = pv.OpenDataFile("../out/fields/{}".format(grid_file))
        jobs = len(start_location)

        print "Looking for Drift Shell at: {} RE".format(start_location[rank])

        # drift_shell = ds.driftShell(PV_dataObject=data, local_times=calcTimes, start_line=(start_location[rank],0,0), K=0, mode='location_k')
        drift_shell = ds.driftShell(local_times=calcTimes, start_line=(start_location[rank],0,0), K=0, mode='analytic_K')

        ideal_RE = ih.mag(drift_shell.field_lines[0].fieldLinePoints_f[0])
        print "[{}] Ideal L*: {}".format(rank, ideal_RE)

        RE_el = 1.0

        test_grids = np.arange(start=50, stop=600, step=50, dtype=np.int)

        # Output to file
        # fileName = "out/error_analysis/{}_lines/l_star_error_analysis_{}_RE_{}_tab.txt".format(lines, grid_file, start_location[rank])
        fileName = "out/error_analysis/{}_lines/l_star_error_analysis_{}_RE_{}_analytic_tab.txt".format(lines, grid_file, start_location[rank])

        f = open(fileName, 'w')
        f.write("# Error Analysis Data for RE = {}\n".format(start_location[rank]))
        f.write("# Author: Joshua Murphy\n")
        f.write("# Time: {}\n".format(dt.datetime.now().time()))
        f.write("# Date: {}\n".format(dt.datetime.now().date()))
        f.write("# Grid File: {}\n\n\n".format(grid_file))

        line = "Grid_Div \t L* \t Error(%) \t Error(diff)\n"
        f.write(line)

        plot_data[grid_file] = col.OrderedDict()

        for size in test_grids:
            L_star = drift_shell.calculate_L_star(RE=RE_el, lat_divs=size, lon_divs=size)
            error = np.abs((ideal_RE - L_star)/ideal_RE * 100)
            diff = ideal_RE - L_star
            plot_data[grid_file][size] = [L_star, error, diff]

            line = "{} \t {} \t {}% \t {}\n".format(size,L_star,error,diff)
            f.write(line)
            print line

        f.close()

        testData = comm.gather(plot_data, root=0)

        if rank == 0:
            data = col.OrderedDict()
            ct = 0
            for r in testData:
                print "Processing Start Location: {}".format(start_location[ct])
                data[start_location[ct]] = col.OrderedDict()
                for file_name in r:
                    for size in r[file_name]:
                        val = r[file_name][size]
                        data[start_location[ct]][size] = val

                ct += 1

            for key in data:
                X,Y = ih.dict_to_x_y(data[key])
                X = [a**2/1000 for a in X]
                Y1 = [a[1] for a in Y]
                Y2 = [a[2] for a in Y]
                Y2 = np.abs(Y2)

                import matplotlib.pylab as plt
                fig1 = plt.figure()
                ax1 = fig1.add_subplot(111)
                ax1.set_ylabel("% Error")
                ax1.set_title("Error Response for L="+str(key)+" vs. Integration Grid Size")

                ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
                ax2.yaxis.tick_right()
                ax2.yaxis.set_label_position("right")
                ax2.set_ylabel("Difference")

                ax1.set_xlabel("Number of Integration Cells ($\\times$ 1000)")

                line1, = ax1.plot(X, Y1, 'bo-', label="% Error")
                line2, = ax2.plot(X, Y2, 'rx-', label="Difference")
                plt.legend(handles=[line1, line2], bbox_to_anchor=(0., 1.02, 1., .102), loc=1, ncol=2, mode="expand", borderaxespad=0.)
                plt.tight_layout()
                # plt.savefig("out/error_analysis/{}_lines/l_star_error_analysis_plot_{}_RE_{}_lines.pdf".format(lines, grid_file, start_location[rank]))
                plt.savefig("out/error_analysis/{}_lines/l_star_error_analysis_analytic_plot_{}_RE_{}_lines.pdf".format(lines, grid_file,
                                                                                               start_location[rank]))



