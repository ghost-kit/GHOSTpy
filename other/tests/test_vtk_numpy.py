# import vtk
import pyevtk.hl as hl

from grids import skelGrid

# from vtk.numpy_interface import dataset_adapter as dsa
#
# import inv_common as inv
#
# reader = vtk.vtkXMLStructuredGridReader()
# reader.SetFileName("test_data/lfm_grid_test.vts")
# reader.Update()
#
# lfm = dsa.WrapDataObject(reader.GetOutput())
# points = np.array(lfm.Points)
#
#
# print points

lfmgrid = skelGrid.Grid(skel_file="test_data/lfm_grid_test.vts")
bdipole = lfmgrid.dipole_field()
bt96 = lfmgrid.t96_field(dyn_pressure=5., dst=-50., imf_by=0., imf_bz=-10., ps=0.0)

hl.gridToVTK("test_data/lfm_dipole_test", lfmgrid.x, lfmgrid.y, lfmgrid.z, pointData={"B": bdipole})
hl.gridToVTK("test_data/lfm_t96_test", lfmgrid.x, lfmgrid.y, lfmgrid.z, pointData={"B": bt96})

print "Done"




