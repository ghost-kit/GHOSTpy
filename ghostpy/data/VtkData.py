import numpy as np
import types as tp
import vtk
from scipy.interpolate import LinearNDInterpolator as gd
from vtk.numpy_interface import dataset_adapter as dsa
import ghostpy.algorithms.FieldTracers as fts

import GpData as gpd
import ghostpy.algorithms.common as algc


class VtkData(gpd.data):
    """
    This class utilizes vtk structured grid .VTK xml files.  To use another type, please
    subclass and change the "vtk_reader()" method
    """
    def __init__(self, filename=None, vector=None, name="vtkData"):
        """
        Initializes a GHOSTpy data object for a vtk data file
        :param filename: Name of the vtk file to read
        :param vector: Name of the vector to be read within the vtk data file
        :param name: common name for the reader. Default is "vtkData"  used to display the trace mode
        """
        assert isinstance(name, tp.StringType), "Must specify a string Name for the data Mode"
        self.name = name
        assert isinstance(filename, tp.StringType), "Filename must be a valid string"
        self.file = filename
        self.reader = self.vtk_reader()
        assert isinstance(vector, tp.StringType), "Vector Name must be a string"
        odata = self.reader.GetOutput()
        self.in_data = dsa.WrapDataObject(odata)
        self.points = self.in_data.GetPoints()
        self.ib = self.find_ib()
        self.array = vector
        self.integrator = fts.VTKrk45(self)
        self.extents = self.reader.GetUpdateExtent()
        self.dims = np.array([0,0,0,0])
        self.dims[0] = self.extents[1]+1
        self.dims[1] = self.extents[3]+1
        self.dims[2] = self.extents[5]+1
        self.dims[3] = 3

        self.grid = self.__build_grid(self.points, self.dims[:-1])

        self.BL = None #self.__bisection_model()


    def __str__(self):
        return "GHOSTpy {} Reader".format(self.name)

    def vtk_reader(self):
        """
        Generates and returns the vtk reader object for the file associated with the data object
        :return: vtkObject for the reader
        """
        reader = vtk.vtkXMLStructuredGridReader()
        reader.SetFileName(self.file)
        reader.Update()
        assert isinstance(reader, vtk.vtkObject)
        return reader

    def get_grid(self):
        return self.grid

    @staticmethod
    def __build_grid(data, dims):
        x, y, z = np.hsplit(data, 3)

        # The grids use Fortran ordering. Need to remember that.
        x = x.reshape(dims, order='F')
        y = y.reshape(dims, order='F')
        z = z.reshape(dims, order='F')

        return x, y, z

    def get_name(self):
        return self.name

    def get_calc_boundary(self):
        return self.ib

    def get_trace_boundary(self):
        return self.get_calc_boundary()

    def get_xyz(self, xyz):

        # FIXME: I actually do have to implement this. Sigh.  I have to use point probes to get the values,
        # FIXME: But I can do them all at once to save time. (this is needed during the resample)

        # 1) define point
        # 2) apply filter
        # 3) return point
        point = vtk.vtkPointSource()
        point.SetCenter(xyz)
        point.SetNumberOfPoints(1)
        point.SetRadius(1e-12)
        point.Update()

        output = self.reader.GetOutput()
        b_field = output.GetPointData().GetArray(self.array)
        output.GetPointData().SetVectors(b_field)
        probe = vtk.vtkProbeFilter()
        probe.SetInputConnection(point.GetOutputPort())
        probe.SetSourceData(output)
        probe.Update()

        pointN = dsa.WrapDataObject(probe.GetOutput())
        pVal = pointN.GetPointData().GetArray(self.array)
        # print ("Value: {}".format(pVal.flatten()))

        return pVal.flatten()

    def get_rlp(self, rlp):
        xyz = algc.sphere_to_cart(r=rlp[0], lam=rlp[1], phi=rlp[2])
        return self.get_xyz(xyz=xyz)

    def get_reader(self):
        return self.reader

    def get_integrator(self):
        return self.integrator

    def __bisection_model(self):
        import scipy.interpolate as interp
        import ghostpy.algorithms.common as algc
        from ghostpy.algorithms import DipoleField as dpf
        ds = np.linspace(start=1.0, stop=55, num=2000)
        B = []
        Ls = []
        for L in ds:
            dipole_L = L
            Ls.append(dipole_L)
            loc = algc.sphere_to_cart(r=dipole_L, lam=0, phi=0)
            B.append(algc.mag(dpf.dipole_field(x=loc[0], y=loc[1], z=loc[2])))

        B = np.array(B, dtype=np.float_)
        return interp.interp1d(B, Ls, kind='cubic')


    def get_bisection_model(self):
        return self.BL


    def find_ib(self):
        gridRE = algc.mag(self.points)
        gridRE = np.around(gridRE, 3) * 1.01
        return np.nanmin(gridRE)