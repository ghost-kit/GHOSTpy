# File:        testFields.py
# -----------------------------------------------------------------------------
# Author:      Joshua J. Murphy
# Project:     PhD dissertation - gridded magnetic field generation
# Purpose:     Generates test fields for T96 and Dipole. Exports the Field Geometry
#              as structured grid VTK files.
# -----------------------------------------------------------------------------

import numpy as np
from pyevtk.hl import gridToVTK

import bfieldgrid

#X = bfieldgrid.bFieldGrid(extents=[20, -60, -20, 20, -20, 20], resolution=(80, 40, 40))
# Y = bfieldgrid.bFieldGrid(extents=[-40, 15, -20, 20, -20, 20], resolution=(56*3, 41*3, 41*3), exclusion_distance=0.2)
# X = bfieldgrid.bFieldGrid(extents=[-40, 15, -20, 20, -20, 20], resolution=(56*10, 41*10, 41*10), exclusion_distance=0.2)
D = bfieldgrid.bFieldGrid(extents=[-15, 15, -15, 15, -15, 15], resolution=(31 * 10, 31 * 10, 31 * 10), exclusion_distance=0.2)
E = bfieldgrid.bFieldGrid(extents=[-15, 15, -15, 15, -15, 15], resolution=(31 * 3, 31 * 3, 31 * 3), exclusion_distance=0.2)


# # Save to disk
# data = X.t96_field(dyn_preassure=5., dst=-50., imf_by=0., imf_bz=-10., dipole_tilt=20.*np.pi/180)
# gridToVTK("./out/fields/t96_dp20_dst50_10_large_grid", X.x, X.y, X.z, pointData={"B": data})
# del data
#
# data = X.t96_field(dyn_preassure=5., dst=-50., imf_by=0., imf_bz=-10., dipole_tilt=340.*np.pi/180)
# gridToVTK("./out/fields/t96_dp340_dst50_10_large_grid", X.x,X.y,X.z, pointData={"B": data})
# del data
#
# data = X.t96_field(dyn_preassure=5., dst=-50., imf_by=0., imf_bz=-10., dipole_tilt=0.0)
# gridToVTK("./out/fields/t96_dp0_dst50_10_large_grid", X.x,X.y,X.z, pointData={"B": data})
# del data

data = D.dipole_field(ps=20*np.pi/180)
gridToVTK("./out/fields/dipole_dp20_10_large_grid", D.x, D.y, D.z, pointData={"B": data})
del data

data = D.dipole_field(ps=0.0)
gridToVTK("./out/fields/dipole_dp0_10_large_grid", D.x, D.y, D.z, pointData={"B": data})
del data

# gridToVTK("./out/fields/t96_dp20_dst50_3_small_grid", Y.x, Y.y, Y.z, pointData={"B": Y.t96_field(dyn_preassure=5., dst=-50., imf_by=0., imf_bz=-10., dipole_tilt=20.*np.pi/180)})
# gridToVTK("./out/fields/t96_dp340_dst50_3_small_grid", Y.x,Y.y,Y.z, pointData={"B": Y.t96_field(dyn_preassure=5., dst=-50., imf_by=0., imf_bz=-10., dipole_tilt=340.*np.pi/180)})
# gridToVTK("./out/fields/t96_dp0_dst50_3_small_grid", Y.x,Y.y,Y.z, pointData={"B": Y.t96_field(dyn_preassure=5., dst=-50., imf_by=0., imf_bz=-10., dipole_tilt=0.0)})
gridToVTK("./out/fields/dipole_dp20_3_small_grid", E.x, E.y, E.z, pointData={"B": E.dipole_field(ps=20*np.pi/180)})
gridToVTK("./out/fields/dipole_dp0_3_small_grid", E.x, E.y, E.z, pointData={"B": E.dipole_field(ps=0.0)})

print "Done."
# print "Origin Exclusion: ", X.origin_exclusion
