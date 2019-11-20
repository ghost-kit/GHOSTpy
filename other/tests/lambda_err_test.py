import matplotlib.pyplot as plt
import numpy as np

from prototype import inv_common as ic


def lambda_intersection_error():
    lvalues = np.linspace(start=1.5, stop=10.5, num=25)

    err_tol = 1e-5

    calc_lambda={}
    anal_lambda={}
    err_lambda={}
    rerr_lambda={}

    for l in lvalues:
        print "Processing L: {}".format(l)

        anal_lambda[l] = ic.analytic_dipole_line_lat_crossing(l,1.0)

        path, value = ic.trace_rk45(x0=[l, 0.0, 0.0], error_tol=err_tol)
        re = []
        for a in path:
            re.append(ic.mag(a))
        re.reverse()
        path.reverse()
        lambda_loc = ic.lin_interp(path, re, 1.0)

        calc_lambda[l] = ic.get_lat(lambda_loc)
        err_lambda[l] = np.abs(anal_lambda[l] - calc_lambda[l])
        rerr_lambda[l] = err_lambda[l]/anal_lambda[l]

        # print "analytic - calc: {} - {} = {}".format(anal_lambda[l], calc_lambda[l], err_lambda[l])
        # print "relative Error: {}".format(rerr_lambda[l])

    Xa, Ya = ic.dict_to_x_y(anal_lambda)
    Xc, Yc = ic.dict_to_x_y(calc_lambda)
    Xe, Ye = ic.dict_to_x_y(err_lambda)
    Xr, Yr = ic.dict_to_x_y(rerr_lambda)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax2 = ax1.twinx()
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")

    line1, = ax1.plot(Xa, Ya, 'ko-', label="Analytic $\lambda$")
    line2, = ax1.plot(Xc, Yc, 'r.-', label="Traced $\lambda$")
    line3, = ax2.semilogy(Xe, Ye, 'b--', label="Absolute Error")

    ax1.set_ylabel("$\lambda$ of Field Line Intersection")
    ax2.set_ylabel("$\lambda$ Error")
    ax1.set_title("$\lambda$ Error vs. $L^{{*}}$ on a dipole\nRK4/5 Error Tol: {}".format(err_tol), y=1.04)
    ax1.set_xlabel("$L^{*} for Dipole$")
    ax2.legend(handles=[line1, line2, line3], loc='upper center', bbox_to_anchor=(0.5, 1.065), ncol=3, fancybox=True, fontsize=9)

    fig1.tight_layout()
    fig1.savefig("out/error_analysis/Lambda_trace_error.pdf")


def e_trace_profile(re= 6.6, plot=False):
    # generate E_trace profile
    Etrace = {}
    Eratio = {}
    lambda_error = np.logspace(start=-14, stop=-1, num=1000, dtype=np.float64)
    lambda_base = ic.analytic_dipole_line_lat_crossing(re, 1.0)

    for err in lambda_error:
        lat = lambda_base + err
        Etrace[err] = np.abs(ic.analytic_dipole_L(lambda_base, r=1.0) - ic.analytic_dipole_L(lat, r=1.0))
        Eratio[err] = Etrace[err]/err

    if plot:
        XeT, YeT = ic.dict_to_x_y(Etrace)
        XeR, YeR = ic.dict_to_x_y(Eratio)

        fig2 = plt.figure()
        ax3 = fig2.add_subplot(111)
        line3, = ax3.loglog(XeT, YeT, 'g-', label="$E_{trace}$", lw=0.5)
        ax3.set_ylabel("$E_{trace}$")
        ax3.set_xlabel("Error in $\lambda$")
        ax3.set_title("Error in $L^*$ due to Error in Field Line Tracing\n($E_{{trace}}$) on a dipole @ $L^*=\ {}\ R_E$".format(re),  y=1.05)
        ax3.grid()

        ax4 = ax3.twinx()
        ax4.yaxis.tick_right()
        ax4.yaxis.set_label_position("right")
        line4,= ax4.plot(XeR, YeR, 'b-', label="Error Effect Ratio", lw=0.5)
        ax4.set_ylabel("Error Effect Ratio")
        ax4.set_ylim([0, 1.5])

        ax4.legend(handles=[line3, line4], loc='upper center', bbox_to_anchor=(0.5, 1.06), ncol=2, fancybox=True)

        # add median label
        mindex = np.argsort(XeR)[len(XeR)//2]
        median_loc = [XeR[mindex], YeR[mindex]]

        ax4.annotate('Median: {:.4f}'.format(YeR[mindex]), xy=median_loc)

        fig2.tight_layout()
        fig2.savefig("out/error_analysis/Etrace/E_trace{}RE.pdf".format(re))
        plt.close(fig2)

    return Etrace, Eratio


def e_trace_profiler(crange=None, num=10, plot=True, subplot=False):
    # default range
    if crange is None:
        crange = [0, 90]

    # get the range
    points = np.linspace(start=crange[0], stop=crange[1], num=num, endpoint=False)
    points2 = ic.analytic_dipole_L(points, 1.0)

    print points2[0]

    profile = {}
    count = 0
    for point in points:
        re = points2[count]
        Et, Er = e_trace_profile(re=re, plot=subplot)
        profile[point] = np.median(Er.values())
        count += 1

    if plot:
        X,Y = ic.dict_to_x_y(profile)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.semilogy(X,Y, 'r.-', label="Median Error Effect Ratio")
        ax.set_xlabel("Base $\lambda$")
        ax.set_ylabel("Error Effect Ratio")
        ax.set_title("Sensitivity of $L^*$ Calculations to Errors in $\lambda$\n in a Dipole")

        # ax2 = ax.twiny()
        # ax2.set_xticklabels(X2)
        # ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
        # ax2.set_xlim(ax.get_xlim())
        # ax2.set_xlabel("Ideal $\lambda$")
        # ax2.xaxis.set_label_position("top")
        ax.grid()

        fig.tight_layout()
        fig.savefig("out/error_analysis/Etrace/E_trace_profile.pdf")
        plt.close(fig)

e_trace_profiler(num=90,subplot=False)
e_trace_profile(re=6.6, plot=True)

print "Done."
