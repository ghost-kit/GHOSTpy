import ghostpy.grids.skelGrid as skel
from pyevtk.hl import gridToVTK
import numpy as np

quad_grid = skel.Grid(skel_file="test_data/WHIQuad.vts")
double_grid = skel.Grid(skel_file="test_data/WHIDouble.vts")
single_grid = skel.Grid(skel_file="test_data/WHISingle.vts")

dipole_quad = quad_grid.dipole_field()
dipole_single = single_grid.dipole_field()
dipole_double = double_grid.dipole_field()

t96_quad = quad_grid.t96_field(dst=-10, imf_by=-5, imf_bz=-10, dyn_pressure=4)
t96_double = double_grid.t96_field(dst=-10, imf_by=-5, imf_bz=-10, dyn_pressure=4)
t96_single = single_grid.t96_field(dst=-10, imf_by=-5, imf_bz=-10, dyn_pressure=4)


print ("Building LFM Quad Grid with dipole")
gridToVTK("test_data/lfm_dipole_test_quad", quad_grid.x, quad_grid.y, quad_grid.z, pointData={"B": dipole_quad})
print ("Building LFM Quad Grid with T96")
gridToVTK("test_data/lfm_t96_test_quad", quad_grid.x, quad_grid.y, quad_grid.z, pointData={"B": t96_quad})

print ("Building LFM Double Grid with dipole")
gridToVTK("test_data/lfm_dipole_test_single", single_grid.x, single_grid.y, single_grid.z, pointData={"B": dipole_single})
print ("Building LFM Double Grid with T96")
gridToVTK("test_data/lfm_t96_test_single", single_grid.x, single_grid.y, single_grid.z, pointData={"B": t96_single})

print ("Building LFM Single Grid with dipole")
gridToVTK("test_data/lfm_dipole_test_double", double_grid.x, double_grid.y, double_grid.z, pointData={"B": dipole_double})
print ("Building LFM Single Grid with T96")
gridToVTK("test_data/lfm_t96_test_double", double_grid.x, double_grid.y, double_grid.z, pointData={"B": t96_double})

print ("DONE")
