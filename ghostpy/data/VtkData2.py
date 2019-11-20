import types as tp
import vtk
from vtk.util import numpy_support as vnp
import numpy as np
import GpData as gpd
import algorithms.common as cm


class VtkData(gpd.data):
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
        self.reader = self.vtk_xml_reader
        assert isinstance(vector, tp.StringType), "Vector Name must be a string"
        self.vector = vector

    @property
    def vtk_xml_reader(self):
        """
        Generates and returns the vtk reader object for the file associated with the data object
        :return: vtkObject for the reader
        """
        reader = vtk.vtkXMLStructuredGridReader()
        reader.SetFileName(self.file)
        reader.Update()
        assert isinstance(reader, vtk.vtkObject)
        return reader

    def get_name(self):
        return self.name

    def get_xyz(self, xyz):
        """
        responds to list of xyz components and returns the values of data at the components
        :param xyz: array of (x,y,z) components to be probed
        :return: array of values coorosponding to (x,y,z) components provided
        """
        if cm.mag(xyz) < self.get_actual_inner_boundary():
            val = np.array([np.NaN, np.NaN, np.NaN])
        else:
            points = self.__get_points_object__([xyz])
            val = self.__get_data_at_points__(points)[0]

        # print (val)
        return val

    def get_rlp(self, rlp):
        """
        responds to list of rlp components and returns the values of data at the components
        :param rlp: array of (r,l,p) components to be probed
        :return: array of values coorosponding to the (r,l,p) components provided
        """
        xyz = cm.sphere_to_cart(r=rlp[0], lam=rlp[1], phi=rlp[2])
        return self.get_xyz(xyz=xyz)

    def get_calc_boundary(self):
        return 3.0

    def get_trace_boundary(self):
        return 2.5

    def get_actual_inner_boundary(self):
        return 2.3

    @staticmethod
    def __get_points_object__(xyz):
        """
        get the vtk point object for the location (x,y,z)
        :param xyz: must be in the form (x1,y1,z1),(x2,y2,z2), ..., (xn, yn, zn)]
        :return: vtk point object
        """
        # TODO: Need to fix this to handle all points, not just the first
        source = vtk.vtkPointSource()
        source.SetCenter(xyz[0])
        source.SetRadius(0)
        source.SetNumberOfPoints(1)
        source.Update()
        return source

    def __get_data_at_points__(self, points):
        """
        extract the values at the given points from the dataset
        :param points: a vtkObject of points
        :return: list of point values
        """
        file_data = self.reader.GetOutput()
        probe = vtk.vtkProbeFilter()
        probe.SetInputConnection(points.GetOutputPort())
        probe.SetSourceData(file_data)
        probe.Update()

        b = vnp.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray(self.vector))
        return b
