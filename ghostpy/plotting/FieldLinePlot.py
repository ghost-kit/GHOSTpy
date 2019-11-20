from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.lines import Line2D
from collections import OrderedDict

from ghostpy.algorithms import common as algc
from ghostpy.Invariants import LShell as ls
import numpy as np
import ghostpy.Invariants.FieldLine as fl
import matplotlib.pylab as plt

class FieldLinePlot(object):

    def __init__(self, _3D_=True):

        self.fig = plt.figure(figsize=(15,15), dpi=150)
        self._3D_ = _3D_

        if self._3D_:
            self.ax = self.fig.add_subplot(111, projection='3d')
        else:
            self.ax = self.fig.add_subplot(111)

        self.ax.set_autoscale_on(False)
        self.set_axes_lim()
        self.set_xlabel("X-Axis")
        self.set_ylabel("Y-Axis")
        if self._3D_:
            self.set_zlabel("Z-Axis")

    def __del__(self):
        plt.close(self.fig)

    def set_axes_lim(self, limits=None):
        if limits is None:
            limits = [[-10, 10], [-10, 10], [-10, 10]]

        self.ax.set_xlim(limits[0])
        self.ax.set_ylim(limits[1])
        if self._3D_:
            self.ax.set_zlim(limits[2])

    def addFieldLine(self, line, color=None, label=""):
        assert isinstance(line, fl.FieldLine)
        trace_n = line.trace_n
        trace_s = line.trace_s

        xn,yn,zn = np.hsplit(trace_n, 3)
        xs,ys,zs = np.hsplit(trace_s, 3)

        if color == None:
            if line.valid:
                color = 'k-'
            else:
                color = 'r-'


        self.ax.plot(xn.flatten(), yn.flatten(), zn.flatten(), color, lw=0.75, label=label)
        self.ax.plot(xs.flatten(), ys.flatten(), zs.flatten(), color, lw=0.75, label=label)

        # plt.show(self.fig)


    def changeCameraAngle(self, theta, phi):
        self.ax.view_init(theta,phi)
        plt.draw()

        pass

    def plot_shell_field_lines(self, lshell, num_lines=10):
        assert isinstance(lshell, ls.LShell)

        if isinstance(lshell.lines, dict):
            print ("We have {} retained lines.".format(lshell.get_number_of_traces()))
            keys = lshell.lines.keys()
            keys = sorted(keys)
            for key in keys:
                if lshell.lines[key][0] is not None:
                    l1 = lshell.lines[key][0]
                    k = lshell.k
                    b = lshell.b
                    assert isinstance(l1, fl.FieldLine)
                    gap = l1.__get_min_b_gap__(b=b, k=k)

                    if gap > l1.error_tol:
                        color = 'r'
                    else:
                        color = 'k'
                    self.addFieldLine(lshell.lines[key][0], color=color, label="LShell Line")


                if lshell.lines[key][1] is not None:
                    l1 = lshell.lines[key][1]
                    k = lshell.k
                    b = lshell.b
                    assert isinstance(l1, fl.FieldLine)
                    gap = l1.__get_min_b_gap__(b=b, k=k)
                    if gap > l1.error_tol:
                        color = 'r'
                    else:
                        color = 'k'
                    self.addFieldLine(lshell.lines[key][1], color=color, label="LShell Line")



        else:
            print ("Generating Drift Trajectory Plot")
            self.add_title("Drift Trajectory for $B_{{mirror}} = {}$, $K = {}$".format(lshell.b, lshell.k))
            drift_shell = lshell.get_drift_trajectory(num_lines)
            for line in drift_shell:
                self.addFieldLine(line=fl.FieldLine(data=lshell.data, start=algc.sphere_to_cart(r=line[0], lam=line[1], phi=line[2])))

        self.ax.scatter(lshell.start[0], lshell.start[1], lshell.start[2], c='r', label="Start Point")

    def plot_drift_boundary(self, lshell):
        if lshell.get_number_of_traces() < 1:
            return None
        assert isinstance(lshell, ls.LShell)
        traj = lshell.get_drift_trajectory_cart(res=100)

        # print (traj)
        x, y, z = np.hsplit(traj, 3)

        x = x.flatten()
        y = y.flatten()
        z = z.flatten()

        if lshell.valid:
            color = 'b-'
        else:
            color = 'r-'

        self.ax.plot(x,y,z, color, alpha=0.85)


    def add_point(self, loc,  marker="o", color="k"):
        x = loc[0]
        y = loc[1]
        z = loc[2]
        self.ax.scatter(x, y, z, marker=marker, color=color)

    def plot_k_b_model(self, line):
        pass

    def add_line(self, start, stop, color='k', lw=2, label=""):
        x = [start[0], stop[0]]
        y = [start[1], stop[1]]
        z = [start[2], stop[2]]

        self.ax.plot(x,y,z, color, lw=lw, label=label)


    def legend(self, loc=1):
        handles, labels = self.ax.get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), loc=loc )

    def show(self):
        # plt.tight_layout()
        plt.show(self.fig)

    def savePDF(self, filename):
        self.fig.savefig(filename)

    def add_title(self, title):
        assert isinstance(title, str)
        self.ax.set_title(title)

    def set_xlabel(self, label):
        assert isinstance(label, str)
        self.ax.set_xlabel(label)

    def set_ylabel(self, label):
        assert isinstance(label, str)
        self.ax.set_ylabel(label)

    def set_zlabel(self, label):
        assert isinstance(label, str)
        if self._3D_:
            self.ax.set_zlabel(label)

    def plot_bk_lines(self, line):
        assert isinstance(line, fl.FieldLine)
        b = line.m_trace_b_mirror
        k = line.K

        x = np.arange(len(b))

        print ("X: {}".format(x))
        print ("b: {}".format(b))
        print ("k: {}".format(k))

        self.ax.semilogy(x, b, label="B")
        # self.ax.semilogy(x, k, label="K")




        self.ax.plot()
