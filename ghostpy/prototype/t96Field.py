# File:        t96Field.py
# -------------------------------------------------------------------------
# Author:      Joshua J. Murphy
# Project:     PhD dissertation dipole grid
# Purpose:     Generates a T-96 magnetic field on a cartesian grid
#                  Saves the grid in vtk structured grid format (.vts)
# -------------------------------------------------------------------------

from pyevtk.hl import gridToVTK

import bfieldgrid as bfg

grid = bfg.bFieldGrid(extents=[-10, 10, -10, 10, -10, 10], resolution=[128, 128, 128])

# Save to disk
gridToVTK("./out/Fields/t96Field128", grid.x, grid.y, grid.z, pointData={"B" : grid.t96_field()})
