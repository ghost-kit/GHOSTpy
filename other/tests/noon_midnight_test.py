# File: noon_midnight_test.py
# Author:  Joshua Murphy
# Project: PhD Dissertation Proof of Concept
# ------------------------------------------------------
# Purpose:
#      Extracts test information from a noon and midnight field line, plots the data, and exports values
#      as ASCII

import matplotlib.pylab as pl
import numpy as np
import paraview.simple as pv
from astropy.io import ascii
from prototype import fieldLine as fl

from ghostpy.prototype import inv_common as ih

pv._DisableFirstRenderCameraReset()

# Load Test Data
t96_128 = pv.OpenDataFile('test_data/lfm_dipole_test.vts')

start_positions = [(-5.5, 0, 0), (-5.25, 0, 0), (-5.00, 0.0, 0.0), (-4.50, 0.0, 0.0), (-4.0, 0.0, 0.0),
                   (-3.5, 0.0, 0.0), (-3.0, 0.0, 0.0), (3.0, 0.0, 0.0), (3.5, 0.0, 0.0),
                   (4.0, 0.0, 0.0), (4.5, 0.0, 0.0), (5.0, 0.0, 0.0), (5.5, 0.0, 0.0),
                   (6.0, 0.0, 0.0), (6.5, 0.0, 0.0),  (6.7, 0.0, 0.0), (6.9, 0.0, 0.0),
                   (7.0, 0.0, 0.0)]
flines = {}
i_integral = {}
dsp = {}
fwd = {}
bkwd = {}
min_sphr = {}
min_sphr_disp = {}
contours = {}
contsphere_b = {}
contsphere_f = {}
rvs = pv.CreateRenderView()
rvs.InteractionMode = '2D'
rvs.CameraPosition = [0, -5, 0]
rvs.CameraFocalPoint = [0, 0, 0]
rvs.CameraViewUp = [0, 0, 1]
rvs.CameraParallelScale = 11
rvs.AxesGrid.Visibility = 1
rvs.AxesGrid.AxesToLabel = 5
rvs.AxesGrid.GridColor = [0.8, 0.8, 0.8]
rvs.AxesGrid.ShowGrid = 1
rvs.AxesGrid.XLabelColor = [0.2, 0.2, 0.2]
rvs.AxesGrid.XTitleColor = [0.2, 0.2, 0.2]
rvs.AxesGrid.ZLabelColor = [0.2, 0.2, 0.2]
rvs.AxesGrid.ZTitleColor = [0.2, 0.2, 0.2]
rvs.OrientationAxesVisibility = 0
rvs.ViewSize = [1920, 1280]
rvs.Background = [1.0, 1.0, 1.0]

# create a sphere to represent Earth
earth = pv.Sphere()
earth.Radius = 1.0
earth.ThetaResolution = 128
earth.PhiResolution = 128
earth_disp = pv.Show(earth, rvs)
earth_disp.DiffuseColor = [0.0, 0.3333333333333333, 0.0]

for start_p in start_positions:
    print "Starting Postion: {}".format(start_p)
    # get field lines and  integrals
    flines[start_p] = fl.fieldLine(t96_128, start=start_p)
    i_integral[start_p] = flines[start_p].get_i_integrals()

    # make field line plots
    dsp[start_p] = pv.Show(t96_128, view=rvs)
    dsp[start_p].CubeAxesVisibility = 0
    dsp[start_p].DiffuseColor = [0.0, 0.0, 0.0]
    dsp[start_p].EdgeColor = [0.0, 0.0, 0.0]
    dsp[start_p].AmbientColor = [0.0, 0.0, 0.0]

    # plot field line
    fwd[start_p] = pv.Show(flines[start_p].fieldLineObj_f,rvs)
    bkwd[start_p] = pv.Show(flines[start_p].fieldLineObj_b,rvs)
    fwd[start_p].DiffuseColor = [0.6666666666666666, 0.0, 0.0]
    bkwd[start_p].DiffuseColor = [0., 0., 0.66]

    # Find B_min
    Bmin_loc, Bmin_val = flines[start_p].get_field_min()

    test_contour_list = []
    test_points_b = []
    test_points_f = []
    test_disp_b = []
    test_disp_f = []
    num_divs = 4
    if not flines[start_p].is_bifurcated:
        for n in range(int(len(flines[start_p].fieldLinePoints_f) / num_divs)):
            test_points_b.append(pv.Sphere())
            test_points_b[-1].Radius = 0.05
            test_points_b[-1].Center = flines[start_p].fieldLinePoints_b[(n * num_divs)]
            test_disp_b.append(pv.Show(test_points_b[-1], rvs))
            test_disp_b[-1].DiffuseColor = [0.1, 0.1, 0.9]

            test_points_f.append(pv.Sphere())
            test_points_f[-1].Radius = 0.05
            test_points_f[-1].Center = flines[start_p].fieldLinePoints_f[(n * num_divs)]
            test_disp_f.append(pv.Show(test_points_f[-1], rvs))
            test_disp_f[-1].DiffuseColor = [0.9, 0.1, 0.1]

            test_contour_list.append(ih.mag(flines[start_p].fieldLineData_f[(n * num_divs)]))

        contours[start_p] = test_contour_list
        contsphere_b[start_p] = test_points_b
        contsphere_f[start_p] = test_points_f

    min_sphr[start_p] = pv.Sphere()
    min_sphr[start_p].Radius = 0.06
    min_sphr[start_p].Center = Bmin_loc
    min_sphr_disp[start_p] = pv.Show(min_sphr[start_p],rvs)
    min_sphr_disp[start_p].DiffuseColor = [0.0, 0.0, 0.0]

Bmag_calc = pv.Calculator(Input=t96_128)
Bmag_calc.ResultArrayName = 'Bmag'
Bmag_calc.Function ='mag(B)'
Bmag_slice = pv.Slice(Input=Bmag_calc)
Bmag_slice.SliceType = 'Plane'
Bmag_slice.SliceOffsetValues = [0.0]
Bmag_slice.SliceType.Origin = [-20.0, 0.0, 0.0]
Bmag_slice.SliceType.Normal = [0.0, 1.0, 0.0]

Bmag_contour = {}

for key in contsphere_b:
    for x in contsphere_f[key]:
        pv.Hide(x, view=rvs)
    for x in contsphere_b[key]:
        pv.Hide(x, view=rvs)

pv.Render(view=rvs)


for key in contours:
    for x in contsphere_f[key]:
        pv.Show(x, view=rvs)
    for x in contsphere_b[key]:
        pv.Show(x, view=rvs)

    Bmag_contour[key] = pv.Contour(Input=Bmag_slice)
    Bmag_contour[key].ContourBy = ['POINTS', 'Bmag']
    Bmag_contour[key].PointMergeMethod = 'Uniform Binning'
    Bmag_contour[key].Isosurfaces = contours[key]
    contour_disp = pv.Show(Bmag_contour[key], rvs)
    contour_disp.DiffuseColor = [0.0, 0.0, 0.0]
    pv.Render(view=rvs)
    pv.WriteImage("verification_" + str(key) + ".png", view=rvs)

    # Hide items not needed for next loop iteration
    pv.Hide(Bmag_contour[key], view=rvs)
    for x in contsphere_f[key]:
        pv.Hide(x, view=rvs)
    for x in contsphere_b[key]:
        pv.Hide(x, view=rvs)
    pv.Render(view=rvs)


# Plot the values of I for each field line
for key in i_integral:
    x, y = ih.dict_to_x_y(i_integral[key])
    fig = pl.figure()
    pl.plot(x,y)
    pl.xlabel("$B_{mirror}$ (nT)")
    pl.ylabel("$I(B_{mirror})$")
    pl.title("Field Line Geometry Integral $I$ with respect to $B_{mirror}$ "
             "\n For Field line that passes through " + str(key))
    fig.savefig("IwrtBm_" + str(key) + ".png")
    pl.close()

# Dump test data to ASCII file.
# 1) create table dictionary
# 2) include X, Y, Z, |B|, Bx, By, Bz, I in dictionary
for key in flines:
    b, i = ih.dict_to_x_y(i_integral[key])
    new_dict = ih.data_to_dict(B=flines[key].get_b_mag_values(), I=np.array(i), loc_f=flines[key].fieldLinePoints_f,
                               B_f=flines[key].fieldLineData_f, loc_b=flines[key].fieldLinePoints_b,
                               B_b=flines[key].fieldLineData_b)

    # 3) write dictionary to ASCII file utilizing astropy.io.ascii
    try:
        ascii.write(new_dict, 'outData' + str(key) + '.tab.txt', format='tab')
    except:
        print "Failed to write: ", key, " -- Invalid Column Length is likely"
        continue
