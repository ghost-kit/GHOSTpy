# File: skelGrid.py
# Purpose: build a grid based on a skeleton from a vts file
# Project: PhD Thesis
# Date: 21 November 2016
import numpy as np
import sys
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from ghostpy.prototype import inv_common as inc
from ghostpy.data import DipoleData as dpd

class Grid:
    def __init__(self, skel_file=None):
        try:
            self.reader = vtk.vtkXMLStructuredGridReader()
            self.reader.SetFileName(skel_file)
            self.reader.Update()
        except:
            print ("Failed to find skeleton file.  Please make sure the file specified is a valid .vts XML file")
            sys.exit(9)

        ## load as a numpy array
        self.skel = dsa.WrapDataObject(self.reader.GetOutput())
        self.extents = self.reader.GetUpdateExtent()
        print ("Extents: {}".format(self.reader.GetUpdateExtent()))
        # self.grid = inc.cm_to_re(np.array(self.skel.Points))
        self.grid = np.array(self.skel.Points)

        self.dims = np.array([0,0,0])
        self.dims[2] = self.extents[1]+1
        self.dims[1] = self.extents[3]+1
        self.dims[0] = self.extents[5]+1

        self.npdims = np.array([0,0,0,0])
        self.npdims[0] = self.dims[0]
        self.npdims[1] = self.dims[1]
        self.npdims[2] = self.dims[2]
        self.npdims[3] = 3

        self.x, self.y, self.z = self.build_grid(self.grid)

        # assert False
        self.dipole = dpd.DipoleData()

    def lfm_data_field(self):
        data = self.skel.PointData.GetArray("B")

        xyz = np.zeros(shape=self.dims)
        x = np.zeros(shape=self.dims)
        y = np.zeros(shape=self.dims)
        z = np.zeros(shape=self.dims)


        count = 0
        for k in range(self.npdims[2]):
            for j in range(self.npdims[1]):
                for i in range(self.npdims[0]):
                    # print ("index: {}".format(count))
                    # print(data[count])
                    xyz[i, j, k] = tuple(data[count])
                    # y[i, j, k] = data[count][1]
                    # z[i, j, k] = data[count][2]
                    count += 1

        return xyz

    def dipole_field(self):

        return inc.dipole_field(self.x, self.y, self.z, ps=0.0)

    def t96_field(self, dyn_pressure, dst, imf_by, imf_bz, ps=0.0):
        return inc.t96_field(x=self.x,
                             y=self.y,
                             z=self.z,
                             p_dyn=dyn_pressure,
                             dst=dst,
                             imf_by=imf_by,
                             imf_bz=imf_bz,
                             ps=ps)

    def build_grid(self, data):
        xyz = np.array(data)
        x,y,z = np.hsplit(xyz,3)
        x = x.flatten()
        y = y.flatten()
        z = z.flatten()

        x = x.reshape(self.dims)
        y = y.reshape(self.dims)
        z = z.reshape(self.dims)

        return x, y, z
