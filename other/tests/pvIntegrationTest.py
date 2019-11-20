# File: pvIntegration.py
# Author:  Joshua Murphy
# Project: PhD Dissertation Proof of Concept
# ------------------------------------------------------
# Purpose:
#      This file provides integration to the ParaView algorithms for this project.  This will be used to access to the
#      VTK data reader as well as the stream-line tracing algorithms.
#      As the project continues, other integrations will be included here.
#      All algorithms will eventually be integrated directly into ParaView and usable from the Paraview GUI.

import matplotlib.pylab as pl
import paraview.simple as pv
from prototype import fieldLine as fl

from ghostpy.prototype import inv_common as ih

# Open the Dipole 128 x 128 x 128 grid file
dipole128 = pv.OpenDataFile("../out/fields/t96_dp0_dst50_small_grid.vts")
fline = fl.fieldLine(dipole128, start=[8.0, 0.0, 0.0])
I = fline.get_i_integrals()

# #################################################
# To active plotting, uncomment everything below #
# #################################################
x, y = ih.dict_to_x_y(I)
fig1 = pl.figure()
pl.plot(x,y)
pl.xlabel("$B_{mirror}$ (nT)")
pl.ylabel("$I(B_{mirror})$")
pl.title("Field Line Geometry Integral $I$ with respect to $B_{mirror}$ \n For Field line that passes through (-40, 0, 0)")
fig1.savefig("IwrtBm.png")

# Create a paraview render view so we can see visual progress.
# pv._DisableFirstRenderCameraReset()
rv128 = pv.CreateRenderView()
rv128.InteractionMode = '2D'
rv128.CameraPosition = [0, -30, 0]
rv128.CameraViewUp = [0, 0, 1]
rv128.CameraParallelScale = 10
rv128.ViewSize = [1280, 1024]
dipoleDisplay = pv.GetDisplayProperties(dipole128, view=rv128)
# pv.Hide(dipole128, view=rv128)

# create a new 'Sphere'
sphere1 = pv.Sphere()
sphere1.Radius = 1.0
sphere1.ThetaResolution = 64
sphere1.PhiResolution = 64
pv.Show(sphere1, rv128)

# get field line

forward = pv.Show(fline.fieldLineObj_f,rv128)
backward = pv.Show(fline.fieldLineObj_b,rv128)
forward.DiffuseColor = [0.6666666666666666, 0.0, 0.0]
backward.DiffuseColor = [0., 0., 0.66]


# Find B_min
Bmin_loc, Bmin_val = fline.get_field_min()

sphere2 = pv.Sphere()
sphere2.Radius = 0.2
sphere2.Center = Bmin_loc
pv.Show(sphere2,rv128)

# render the views when needed
pv.Render()

# Keep displays up until user hits enter.
raw_input("Press Enter to Continue...")


fline.recompute_field_line(new_start_location=[-5.0, 0.0, 0.0])
Bmin_loc, Bmin_val = fline.get_field_min()
sphere2.Center = Bmin_loc

pv.Render()


# Keep displays up until user hits enter.
raw_input("Press Enter to Continue...")
