import collections as col
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

from prototype import inv_common as ih

lats = np.linspace(start=0, stop=90, num=45, endpoint=False)
data = col.OrderedDict()
data2 = col.OrderedDict()
analytic = col.OrderedDict()
error = col.OrderedDict()
error2 = col.OrderedDict()
absErr = col.OrderedDict()
absErr2 = col.OrderedDict()

divs = np.logspace(start=0.5, stop=3.125, num=25, dtype=int)

maxErr = col.OrderedDict()
maxErr2 = col.OrderedDict()
maxAbs = col.OrderedDict()
maxAbs2 = col.OrderedDict()

for div in divs:
    print "Processing with Cell Count: {}".format(div**2)
    maxErr[div] = []
    maxErr2[div] = []
    maxAbs[div] = []
    maxAbs2[div] = []

    for lat in lats:
        grid = ih.build_spherical_cap_grid(RE=1, divs=div, base_lat=lat, grid_type='sphere')
        data[lat] = ih.integrate_flux_cap_area(sphere_grid=grid)
        data2[lat] = ih.integrate_flux_cap_area2(sphere_grid=grid)
        analytic[lat] = ih.analytic_sphere_cap_area(lat=lat, r=1)
        error[lat] = np.abs((analytic[lat] - data[lat])/analytic[lat])
        error2[lat] = np.abs((analytic[lat] - data2[lat])/analytic[lat])
        absErr[lat] = np.abs(analytic[lat] - data[lat])
        absErr2[lat] = np.abs(analytic[lat] - data2[lat])
        maxErr[div].append(np.max(error[lat]))
        maxAbs[div].append(np.max(absErr[lat]))
        maxErr2[div].append(np.max(error2[lat]))
        maxAbs2[div].append(np.max(absErr2[lat]))


    maxErr[div] = np.max(maxErr[div])
    maxErr2[div] = np.max(maxErr2[div])
    maxAbs[div]= np.max(maxAbs[div])
    maxAbs2[div] = np.max(maxAbs2[div])

    Xabs, Yabs = ih.dict_to_x_y(absErr)
    Xabs2, Yabs2 = ih.dict_to_x_y(absErr2)
    Xfrac, Yfrac = ih.dict_to_x_y(error)
    Xfrac2, Yfrac2 = ih.dict_to_x_y(error2)
    Xcalc, Ycalc = ih.dict_to_x_y(data)
    Xcalc2, Ycalc2 = ih.dict_to_x_y(data2)
    Xanal, Yanal = ih.dict_to_x_y(analytic)

    fig, ax = plt.subplots()

    axes = [ax, ax.twinx(), ax.twinx()]

    fig.subplots_adjust(right=0.750)
    fig.set_size_inches(18.5, 10.5, forward=True)

    axes[-1].spines['right'].set_position(('axes', 1.2))

    axes[-1].set_frame_on(True)
    axes[-1].patch.set_visible(False)

    line1, = axes[0].plot(Xcalc, Ycalc, 'k*-', label="Integration Method 1")
    line2, = axes[0].plot(Xanal, Yanal, 'b.-', label="Analytic Area")
    line3, = axes[1].plot(Xfrac, Yfrac, 'r--', label="Relative Error")
    line4, = axes[2].plot(Xabs, Yabs, 'k-.', label="Absolute Error")

    axes[0].set_title("Spherical Cap Integration (Method 1) Error Profile for ({} x {}) Cells".format(div, div), y=1.04)
    axes[0].set_ylabel("Area")
    axes[1].set_ylabel("Relative Error", color='r')
    axes[2].set_ylabel("Absolute Error")
    axes[0].set_xlabel("Starting $\lambda$")

    axes[1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    axes[2].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    for tl in axes[1].get_yticklabels():
        tl.set_color('r')

    for tl in axes[0].get_yticklabels():
        tl.set_color('b')

    plt.legend(handles=[line1, line2, line3, line4], loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True)
    fig.savefig("out/Spherical_cap_integration_{}.pdf".format(div))
    plt.close(fig)
    fig2, ax2 = plt.subplots()

    axes2 = [ax2, ax2.twinx(), ax2.twinx()]

    fig2.subplots_adjust(right=0.750)
    fig2.set_size_inches(18.5, 10.5, forward=True)

    axes2[-1].spines['right'].set_position(('axes', 1.2))

    axes2[-1].set_frame_on(True)
    axes2[-1].patch.set_visible(False)

    line5, = axes2[0].plot(Xcalc2, Ycalc2, 'k*-', label="Integration Method 2")
    line6, = axes2[0].plot(Xanal, Yanal, 'b.-', label="Analytic Area")
    line7, = axes2[1].plot(Xfrac2, Yfrac2, 'r--', label="Relative Error")
    line8, = axes2[2].plot(Xabs2, Yabs2, 'k-.', label="Absolute Error")

    axes2[0].set_title("Spherical Cap Integration (Method 2) Error Profile for ({} x {}) Cells".format(div, div), y=1.04)
    axes2[0].set_ylabel("Area")
    axes2[1].set_ylabel("Relative Error", color='r')
    axes2[2].set_ylabel("Absolute Error")
    axes2[0].set_xlabel("Starting $\lambda$")

    axes2[1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    axes2[2].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    for tl in axes2[1].get_yticklabels():
        tl.set_color('r')

    for tl in axes2[0].get_yticklabels():
        tl.set_color('b')

    plt.legend(handles=[line5, line6, line7, line8], loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True)
    fig2.savefig("out/Spherical_cap_integration2_{}.pdf".format(div))
    plt.close(fig2)

Xe,Ye = ih.dict_to_x_y(maxErr)
Xe2, Ye2 = ih.dict_to_x_y(maxErr2)
Xa,Ya = ih.dict_to_x_y(maxAbs)
Xa2, Ya2 = ih.dict_to_x_y(maxAbs2)

Xe *= np.array(Xe)
Xe2 *= np.array(Xe2)
Xa *= np.array(Xa)
Xa2 *= np.array(Xa2)

fig2 = plt.figure()
fig2.set_size_inches(18.5, 10.5, forward=True)
ax2 = fig2.add_subplot(111)
lineA, = ax2.loglog(Xe, Ye, 'k*-',label="Max Relative Error (Method 1)")
lineB, = ax2.loglog(Xe2, Ye2,'k*-.', label="Max Relative Error (Method 2)")
lineC, = ax2.loglog(Xa, Ya, 'r.-', label="Max Absolute Error (Method 1)")
lineD, = ax2.loglog(Xa2, Ya2, 'r.-.', label="Max Absolute Error (Method 2)")
ax2.set_title("Max Absolute and Relative Error")
plt.legend(handles=[lineA, lineB, lineC, lineD], loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True)
ax2.set_title("Max Absolute and Relative Error")
ax2.set_ylabel("Error", color='k')
ax2.set_xlabel("Number of Grid Cells")
ax2.set_title("Convergence of Two Spherical Cap Integration Methods", y=1.04)

fig2.savefig("out/convergence.pdf")
plt.close(fig2)
# plt.show()



