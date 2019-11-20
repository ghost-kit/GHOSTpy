# File:        dipole.py
# -------------------------------------------------------------------------
# Author:      Joshua J. Murphy
# Project:     PhD dissertation dipole grid
# Purpose:     Generates a dipole magnetic field on a cartesian grid
#                  Saves the grid in vtk structured grid format (.vts)
# -------------------------------------------------------------------------

from pyevtk.hl import gridToVTK

import bfieldgrid as bfg

grid = bfg.bFieldGrid(extents=[20, -40, -10, 10, -10, 10], resolution=[128, 128, 128])

print "Saving Dipole Field to file...."

# Save to disk
gridToVTK("./out/Fields/dipole128", grid.x, grid.y, grid.z, pointData={"B" : grid.dipole_field()})
