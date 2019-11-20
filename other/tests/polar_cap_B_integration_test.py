import collections as col
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

from prototype import inv_common as ih

reslat = 12
resdivs = 20
stopdivpow = 3.0

lats = np.linspace(start=45, stop=75, num=reslat, endpoint=False)
data = col.OrderedDict()
data2 = col.OrderedDict()
analytic = col.OrderedDict()
error = col.OrderedDict()
error2 = col.OrderedDict()
absErr = col.OrderedDict()
absErr2 = col.OrderedDict()

divs = np.logspace(start=0.01, stop=stopdivpow, num=resdivs, dtype=int)
divs = np.unique(divs)

medianErr = col.OrderedDict()
medianErr2 = col.OrderedDict()
median = col.OrderedDict()
medianAbs2 = col.OrderedDict()

funs = [ih.cell_avg_dipole, ih.cell_cntr_dipole]

for fun in funs:
    medianErr[fun.__name__] = {}
    medianErr2[fun.__name__] = {}
    median[fun.__name__] = {}
    medianAbs2[fun.__name__] = {}
    data[fun.__name__] = {}
    data2[fun.__name__] = {}
    analytic[fun.__name__] = {}
    error[fun.__name__] = {}
    error2[fun.__name__] = {}
    absErr[fun.__name__] = {}
    absErr2[fun.__name__] = {}

    for div in divs:
        print "Processing with Cell Count: {} ({} x {})".format(div**2,div, div)
        medianErr[fun.__name__][div] = 0
        medianErr2[fun.__name__][div] = 0
        median[fun.__name__][div] = 0
        medianAbs2[fun.__name__][div] = 0

        for lat in lats:
            # lets get the integration done
            grid = ih.build_spherical_cap_grid(RE=1, divs=div, base_lat=lat, grid_type='sphere')
            data[fun.__name__][lat] = ih.L_star_from_flux(flux=ih.integrate_flux_cap_area(sphere_grid=grid, pvalue_fun=fun))
            data2[fun.__name__][lat] = ih.L_star_from_flux(flux=ih.integrate_flux_cap_area2(sphere_grid=grid, pvalue_fun=fun))
            analytic[fun.__name__][lat] = ih.analytic_dipole_L(lat, 1)
            error[fun.__name__][lat] = np.abs((analytic[fun.__name__][lat] - data[fun.__name__][lat])/analytic[fun.__name__][lat])
            error2[fun.__name__][lat] = np.abs((analytic[fun.__name__][lat] - data2[fun.__name__][lat])/analytic[fun.__name__][lat])
            absErr[fun.__name__][lat] = np.abs(analytic[fun.__name__][lat] - data[fun.__name__][lat])
            absErr2[fun.__name__][lat] = np.abs(analytic[fun.__name__][lat] - data2[fun.__name__][lat])

        medianErr[fun.__name__][div] = np.median(error[fun.__name__].values())
        median[fun.__name__][div] = np.median(absErr[fun.__name__].values())
        medianErr2[fun.__name__][div] = np.median(error2[fun.__name__].values())
        medianAbs2[fun.__name__][div] = np.median(absErr2[fun.__name__].values())

        # medianErr[fun.__name__][div] = np.median(medianErr[fun.__name__][div])
        # medianErr2[fun.__name__][div] = np.median(medianErr2[fun.__name__][div])
        # median[fun.__name__][div] = np.median(median[fun.__name__][div])
        # medianAbs2[fun.__name__][div] = np.median(medianAbs2[fun.__name__][div])

        Xabs, Yabs = ih.dict_to_x_y(absErr[fun.__name__])
        Xabs2, Yabs2 = ih.dict_to_x_y(absErr2[fun.__name__])
        Xfrac, Yfrac = ih.dict_to_x_y(error[fun.__name__])
        Xfrac2, Yfrac2 = ih.dict_to_x_y(error2[fun.__name__])
        Xcalc, Ycalc = ih.dict_to_x_y(data[fun.__name__])
        Xcalc2, Ycalc2 = ih.dict_to_x_y(data2[fun.__name__])
        Xanal, Yanal = ih.dict_to_x_y(analytic[fun.__name__])

        fig, ax = plt.subplots()

        axes = [ax, ax.twinx(), ax.twinx()]

        fig.subplots_adjust(right=0.750)
        fig.set_size_inches(18.5, 10.5, forward=True)

        axes[-1].spines['right'].set_position(('axes', 1.2))

        axes[-1].set_frame_on(True)
        axes[-1].patch.set_visible(False)

        line1, = axes[0].plot(Xcalc, Ycalc, 'bo-', label="$L^*$ (Method 1)")
        line2, = axes[0].plot(Xanal, Yanal, 'g-', label="Analytic $L^*$")
        line3, = axes[1].plot(Xfrac, Yfrac, 'r-.', label="Relative Error")
        line4, = axes[2].plot(Xabs, Yabs, 'k--', label="Absolute Error")

        axes[0].set_title("$L^*$ (Integration Method 1) Error Profile for ({} x {}) Cells\n(Value Function: {})".format(div, div, fun.__name__), y=1.04)
        axes[0].set_ylabel("$L^*$", color='b')
        axes[1].set_ylabel("Relative Error", color='r')
        axes[2].set_ylabel("Absolute Error")
        axes[0].set_xlabel("Starting $\lambda$")

        axes[1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.4e'))
        axes[2].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.4e'))

        plt.legend(handles=[line1, line2, line3, line4], loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True)
        fig.savefig("out/lstar_m1_error_profile_{}x{}_fun_{}.pdf".format(div,div, fun.__name__))
        plt.close(fig)

        fig2, ax2 = plt.subplots()

        axes2 = [ax2, ax2.twinx(), ax2.twinx()]

        fig2.subplots_adjust(right=0.750)
        fig2.set_size_inches(18.5, 10.5, forward=True)

        axes2[-1].spines['right'].set_position(('axes', 1.2))

        axes2[-1].set_frame_on(True)
        axes2[-1].patch.set_visible(False)

        line5, = axes2[0].plot(Xcalc2, Ycalc2, 'bo-', label="$L^*$ (Method 2)")
        line6, = axes2[0].plot(Xanal, Yanal, 'g-', label="Analytic $L^*$")
        line7, = axes2[1].plot(Xfrac2, Yfrac2, 'r-.', label="Relative Error")
        line8, = axes2[2].plot(Xabs2, Yabs2, 'k--', label="Absolute Error")

        axes2[0].set_title("$L^*$ (Integration Method 2) Error Profile for ({} x {}) Cells\n(Value Function: {})".format(div, div, fun.__name__), y=1.04)
        axes2[0].set_ylabel("$L^*$", color='b')
        axes2[1].set_ylabel("Relative Error", color='r')
        axes2[2].set_ylabel("Absolute Error")
        axes2[0].set_xlabel("Starting $\lambda$")

        axes2[1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.4e'))
        axes2[2].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.4e'))

        plt.legend(handles=[line5, line6, line7, line8], loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True)
        fig2.savefig("out/lstar_m2_error_profile_{}x{}_fun_{}.pdf".format(div,div, fun.__name__))
        plt.close(fig2)

fig2 = plt.figure()
fig2.set_size_inches(18.5, 10.5, forward=True)
ax2 = fig2.add_subplot(111)

fun_names = medianErr.keys()

colors = [['k', 'r'], ['b', 'g'], ['y', 'brown']]

count = 0
for fun_name in fun_names:
    Xe,Ye = ih.dict_to_x_y(medianErr[fun_name])
    Xe2, Ye2 = ih.dict_to_x_y(medianErr2[fun_name])
    Xa,Ya = ih.dict_to_x_y(median[fun_name])
    Xa2, Ya2 = ih.dict_to_x_y(medianAbs2[fun_name])

    Xe *= np.array(Xe)
    Xe2 *= np.array(Xe2)
    Xa *= np.array(Xa)
    Xa2 *= np.array(Xa2)

    pattern1 = '{}*-'.format(colors[count][0])
    pattern2 = '{}.-'.format(colors[count][0])
    pattern3 = '{}*-.'.format(colors[count][1])
    pattern4 = '{}.-.'.format(colors[count][1])

    lineC, = ax2.loglog(Xa, Ya, pattern1, label="Median Error - Area Meth: Sqr - Val Fun: {}".format(fun_name))
    lineD, = ax2.loglog(Xa2, Ya2, pattern4, label="Median Error - Area Meth: Tri - Val Fun: {}".format(fun_name))
    ax2.set_ylabel("Error", color='k')
    ax2.set_xlabel("Number of Grid Cells")
    ax2.set_title("Convergence of $L^*$ (Two Different Methods)", y=1.04)

    count += 1
    # plt.show()

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True)
plt.ylim([1e-8, 1e2])
fig2.savefig("out/convergence_lstar.pdf")
plt.close(fig2)



