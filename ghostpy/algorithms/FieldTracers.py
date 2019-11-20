# Field Line Tracers
import numpy as np
from collections import OrderedDict
from scipy.integrate import ode

import ghostpy.algorithms.common as algc
from ghostpy.data import GpData as gpd
import vtk
from vtk.numpy_interface import dataset_adapter as dsa


class GpIntegrator(object):
    def __init__(self, data=None):
        self.data = data
        assert isinstance(data, gpd.data), "INTEGRATOR DATA SOURCE ERROR: Incorrect Type"
        self.trace_boundary = self.data.get_trace_boundary()

    def integrate(self, x0=None, max_steps=None, error_tol=1e-6, direct='f'):
        pass


class SPrk45(GpIntegrator):
    """
    This class is for integrating analytic fields that have points defined at all points.
    It it can be used on gridded data, but the VTKrk45 class performs much better for gridded data sets
    """

    def __init__(self, data=None):
        GpIntegrator.__init__(self, data)
        self.stop = False

    def sp_rk45(self,  x0=None, max_steps=10000, error_tol=1e-6, direct="f"):
        self.stop = False

        trace = OrderedDict()
        value = OrderedDict()

        assert x0 is not None

        dt = 1e-4
        solver = ode(self.int_fun)
        solver.set_f_params(direct)
        solver.set_integrator('dopri5', rtol=0, atol=error_tol, max_step=0.01, nsteps=max_steps)
        solver.set_solout(self.solout)
        solver.set_initial_value(y=x0)

        value[0] = self.int_fun(0, x0)
        trace[0] = x0

        while solver.successful() and solver.t < max_steps and not self.stop:
            solver.integrate(t=solver.t + dt)
            value[solver.t] = self.int_fun(0, solver.y)
            trace[solver.t] = solver.y.copy()

        # print ("Returning...")
        return np.array(trace.values()), np.array(value.values())

    def solout(self, t, y):
        if algc.mag(y) <= self.trace_boundary:
            self.stop = True
            # print ("location: {}".format(y))
            return -1

    def int_fun(self, t, y, direction='f'):
        val = self.data.get_xyz(y)
        if direction == 'b' and val is not None:
            val *= -1
        return val

    def integrate(self, x0=None, max_steps=10000, error_tol=1e-6, direct='f'):
        return self.rk45(val_fun=self.valfun_dipole, x0=x0, max_steps=max_steps, error_tol=error_tol, direct=direct)
        # return self.sp_rk45(x0=x0, max_steps=max_steps, error_tol=error_tol, direct=direct)


    def valfun_dipole(self, xyz):
        return self.data.get_xyz(xyz=xyz)


    def rk45(self,inner_boundary=0.5, val_fun=None, h=1e-6, x0=None, max_steps=2000, error_tol=1e-3, direct="f"):
        """
        Implementation of the Runge-Kutta 4/5 algorithm
        :param inner_boundary: Where to stop when approaching origin
        :param val_fun: function that can provide a value at (xyz)
        :param h: starting step size
        :param x0: starting point
        :param max_steps: maximum number of steps to compute
        :param error_tol: Error tolerance for calculating step size
        :param direct: Direction of integration ('f' for forward, 'b' for backward)
        :return:
        """

        # print ("Getting Trace:")
        if x0 is None:
            x0 = [6.0, 0.0, 0.0]
        if direct == 'b':
            d = -1
        else:
            d = 1
        mag_x0 = algc.mag(x0)
        x = x0
        RE = mag_x0
        # print ("RE: {}".format(RE))
        dpv = val_fun(x)
        path = [tuple(x)]
        value = [tuple(dpv)]
        steps = 0
        # print ("about to start loop")
        # print ("Inner Boundary: {}".format(inner_boundary))
        while RE > inner_boundary:
            # print ("Inn Loop")
            s = 0
            # print ("S: {}".format(s))
            # while not np.isclose(s, 1.0, atol=0.0001, rtol=0.00) and not np.isnan(s):
            k1 = h * (val_fun(x) * d)
            k2 = h * (val_fun(x + 1. / 4. * k1) * d)
            k3 = h * (val_fun(x + 3. / 32. * k1 + 9. / 32. * k2) * d)
            k4 = h * (val_fun(x + 1932. / 2197. * k1 - 7200. / 2197. * k2 + 7296. / 2197. * k3) * d)
            k5 = h * (val_fun(x + 439. / 216. * k1 - 8. * k2 + 3680. / 513. * k3 - 845. / 4104. * k4) * d)
            k6 = h * (val_fun(x - 8. / 27. * k1 + 2. * k2 - 3544. / 2565. * k3 + 1859. / 4104. * k4 - 11. / 40. * k5) * d)

            x1 = x + 25. / 216. * k1 + 1408. / 2565. * k3 + 2197. / 4104. * k4 - 1. / 5. * k5
            x2 = x + 16. / 135. * k1 + 6656. / 12825. * k3 + 28561. / 56430. * k4 - 9. / 50. * k5 + 2. / 55. * k6

            s = 0.84 * (error_tol * h / (algc.mag(x2 - x1))) ** 0.25

            if np.isnan(s) or np.isinf(s):
                s = 1.0
            h *= s
            x = x1
            RE = algc.mag(x)
            dpv = val_fun(x)
            steps += 1
            if steps > max_steps:
                print("Truncating on number of steps.\nConsider Increasing number of steps.")
                break
            path.append(tuple(x1))
            value.append(tuple(dpv))

        # print ("Path: {}".format(path))

        return np.array(path), np.array(value)






class VTKrk45(GpIntegrator):

    def __init__(self, data=None):
        GpIntegrator.__init__(self, data)
        vtk.vtkObject.GlobalWarningDisplayOff()

    def vtk_rk45(self, x0=None, max_steps=10000, error_tol=1e-6, direct="f"):

        reader = self.data.get_reader()
        assert isinstance(reader, vtk.vtkObject)
        reader.Update()

        streamer1 = self.getStreamer(reader, x0=x0, direct=direct, error_tol=error_tol, maxSteps=max_steps)

        lineN = dsa.WrapDataObject(streamer1.GetOutput())
        pointsN = lineN.GetPoints()
        dataN = lineN.GetPointData().GetArray(self.data.array)

        return pointsN, dataN

    @staticmethod
    def getStreamer(vtkReader, array='B', x0=None, maxSteps=50000, error_tol=1e-8, direct='f'):

        if direct == 'f':
            direction = 0
        elif direct == 'b':
            direction = 1
        else:
            print("Unknown Direction Specified. Please use 'f' for forward (norther hemisphere) "
                  "and 'b' for backward (southern hemisphere)")
            assert False
        # print ("X0: {}".format(x0))
        # There seems to be a problem with the interpolator getting some point exactly on the x axis.
        #   This small oscillation seems to fix the issue.

        x0 += [0,1e-12, -1e-12]

        output = vtkReader.GetOutput()
        b_field = output.GetPointData().GetArray(array)
        output.GetPointData().SetVectors(b_field)
        rk45 = vtk.vtkRungeKutta45()
        streamer = vtk.vtkStreamTracer()
        streamer.SetInputConnection(vtkReader.GetOutputPort())
        streamer.SetInputData(output)
        streamer.SetStartPosition(x0)
        streamer.SetMaximumPropagation(maxSteps)
        streamer.SetIntegrationStepUnit(2)
        streamer.SetMinimumIntegrationStep(0.01)
        streamer.SetMaximumIntegrationStep(2.0)
        streamer.SetInitialIntegrationStep(0.5)
        streamer.SetIntegrationDirection(direction)
        streamer.SetIntegrator(rk45)
        streamer.SetMaximumError(error_tol/10)
        streamer.Update()
        # print("Streamer Complete")
        return streamer

    def integrate(self, x0=None, max_steps=50000, error_tol=1e-6, direct='f'):
        return self.vtk_rk45(x0=x0, max_steps=max_steps, error_tol=error_tol, direct=direct)


