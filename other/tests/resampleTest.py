# File: resampleTest.py
# Author: Joshua Murphy
# Project: PhD Dissertation Proof of Concept
# -----------------------------------------------------------------------------------
# Purpose:
#       This script test the resampling of the fieldLine routine to ensure that it
#       is resampling on the field line.
#
# -----------------------------------------------------------------------------------

import paraview.simple as pv
from prototype import fieldLine as fl

from ghostpy.prototype import inv_common as ih

pv._DisableFirstRenderCameraReset()

# Open T96 grid with dipole tilt
t96_128_dp = pv.OpenDataFile("test_data/lfm_dipole_test.vts")

# Create Paraview render View
rv = pv.CreateRenderView()
rv.InteractionMode = '2D'
rv.CameraPosition = [-20, -30, 0]
#rv.CameraPosition = [-20, 0, 30]

rv.CameraFocalPoint = [-20, 0, 0]
#rv.CameraFocalPoint = [0, 0, 0]

rv.CameraViewUp = [0, 0, 1]
#rv.CameraViewUp = [0, 1, 0]

rv.CameraParallelScale = 23
rv.AxesGrid.Visibility = 1
rv.AxesGrid.AxesToLabel = 5
rv.AxesGrid.GridColor = [0.8, 0.8, 0.8]
rv.AxesGrid.ShowGrid = 1
rv.AxesGrid.XLabelColor = [0.2, 0.2, 0.2]
rv.AxesGrid.XTitleColor = [0.2, 0.2, 0.2]
rv.AxesGrid.ZLabelColor = [0.2, 0.2, 0.2]
rv.AxesGrid.ZTitleColor = [0.2, 0.2, 0.2]
rv.OrientationAxesVisibility = 0

# rv.ViewSize = [1280*3, 600*3]
rv.ViewSize = [1280, 600]

td = pv.Show(t96_128_dp, view=rv)
td.DiffuseColor = [0.0, 0.0, 0.0]
td.EdgeColor = [0.0, 0.0, 0.0]
td.AmbientColor = [0.0, 0.0, 0.0]


# create a sphere to represent Earth
earth = pv.Sphere()
earth.Radius = 1.0
earth.ThetaResolution = 128
earth.PhiResolution = 128
earth_disp = pv.Show(earth, rv)
earth_disp.DiffuseColor = [0.0, 0.3333333333333333, 0.0]

# create test field line
fline = fl.fieldLine(t96_128_dp, start=[-10.0,0.0,0.0])
forward = pv.Show(fline.fieldLineObj_f,rv)
backward = pv.Show(fline.fieldLineObj_b,rv)
forward.DiffuseColor = [0.5, 0.0, 0.0]
backward.DiffuseColor = [0., 0.0, 0.5]

# Mark New Center (Should be B_min)
new_center = fline.fieldLinePoints_f[0]
print "New Center: ", new_center

Bmag_calc = pv.Calculator(Input=t96_128_dp)
Bmag_calc.ResultArrayName = 'Bmag'
Bmag_calc.Function ='mag(B)'
Bmag_slice = pv.Slice(Input=Bmag_calc)
Bmag_slice.SliceType = 'Plane'
Bmag_slice.SliceOffsetValues = [0.0]
Bmag_slice.SliceType.Origin = [-20.0, 0.0, 0.0]
Bmag_slice.SliceType.Normal = [0.0, 1.0, 0.0]
Bmag_contour = pv.Contour(Input = Bmag_slice)
Bmag_contour.ContourBy = ['POINTS', 'Bmag']
Bmag_contour.PointMergeMethod = 'Uniform Binning'

# point = fline.get_location_for_RE(1.0)
#
# endPoint = pv.Sphere()
# endPoint.Radius = 0.05
# endPoint.Center = point[0]
# pv.Show(endPoint,rv)

test_contour_list = []
test_points_b = []
test_points_f = []
test_disp_b = []
test_disp_f = []
for n in range(int(len(fline.fieldLinePoints_f)/10)):
    test_points_b.append(pv.Sphere())
    test_points_b[-1].Radius = 0.20
    test_points_b[-1].Center = fline.fieldLinePoints_b[(n * 10)]
    test_disp_b.append(pv.Show(test_points_b[-1], rv))
    test_disp_b[-1].DiffuseColor = [0.1, 0.1, 0.9]

    test_points_f.append(pv.Sphere())
    test_points_f[-1].Radius = 0.20
    test_points_f[-1].Center = fline.fieldLinePoints_f[(n * 10)]
    test_disp_f.append(pv.Show(test_points_f[-1], rv))
    test_disp_f[-1].DiffuseColor = [0.9, 0.1, 0.1]

    test_contour_list.append(ih.mag(fline.fieldLineData_f[(n*10)]))

# min point
test_points_m = pv.Sphere()
test_points_m.Radius = 0.25
test_points_m.Center = fline.fieldLinePoints_b[0]
test_disp_m = pv.Show(test_points_m, rv)
test_disp_m.DiffuseColor = [0, 0, 0]


Bmag_contour.Isosurfaces = test_contour_list
contour_disp = pv.Show(Bmag_contour, rv)
contour_disp.DiffuseColor = [0.0, 0.0, 0.0]

rv.Background = [1.0, 1.0, 1.0]

pv.Render(rv)

pv.WriteImage("test.png", view=rv)



# Keep displays up until user hits enter.
# raw_input("Press Enter to Continue...")