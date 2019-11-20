import numpy as np
import paraview.simple as pv

from prototype import driftShell as ds

pv._DisableFirstRenderCameraReset()

checkTimes = np.linspace(0, 24, num=5, endpoint=False)
calcTimes = [12.0, 6.0, 18.0, 24.0]

t96 = pv.OpenDataFile("../out/fields/t96_dp0_dst50_small_grid.vts")
t96_K0_ds = ds.driftShell(PV_dataObject=t96, local_times=calcTimes, b_mirror=200.0, K=100, mode='K')
t96_K0_ds2 = ds.driftShell(PV_dataObject=t96, local_times=checkTimes, b_mirror=200.0, K=100, mode='K')

RE_el = 1.0

phi_location = t96_K0_ds.get_phi_locations(RE=RE_el)
print "Locations: ", phi_location

interp_points = np.linspace(0,24,50, endpoint=False)

for x in interp_points:
    if x not in phi_location.keys():
        phi_location[x] = t96_K0_ds.get_new_phi_for_localtime(x, RE=RE_el)


# sys.exit()

rvs4 = pv.CreateRenderView()
rvs4.InteractionMode = '2D'
rvs4.CameraPosition = [0, 55, 0]
rvs4.CameraFocalPoint = [0, 0, 0]
rvs4.CameraViewUp = [0, 1, 0]
rvs4.AxesGrid.Visibility = 1
#rvs4.AxesGrid.AxesToLabel = 5
rvs4.AxesGrid.GridColor = [0.8, 0.8, 0.8]
rvs4.AxesGrid.ShowGrid = 1
rvs4.AxesGrid.XLabelColor = [0.2, 0.2, 0.2]
rvs4.AxesGrid.XTitleColor = [0.2, 0.2, 0.2]
rvs4.AxesGrid.ZLabelColor = [0.2, 0.2, 0.2]
rvs4.AxesGrid.ZTitleColor = [0.2, 0.2, 0.2]
rvs4.AxesGrid.YLabelColor = [0.2, 0.2, 0.2]
rvs4.AxesGrid.YTitleColor = [0.2, 0.2, 0.2]
rvs4.OrientationAxesVisibility = 0
rvs4.ViewSize = [1920, 1280]
rvs4.Background = [1.0, 1.0, 1.0]

# create a sphere to represent Earth
earth = pv.Sphere()
earth.Radius = 1.0
earth.ThetaResolution = 128
earth.PhiResolution = 128

earth_disp4 = pv.Show(earth, rvs4)
earth_disp4.DiffuseColor = [0.51, 0.51, 0.51]
# pv.Hide(earth,rvs4)

for key in t96_K0_ds2.field_lines:
    if t96_K0_ds2.field_lines[key] is None:
        continue
    forward = pv.Show(t96_K0_ds2.field_lines[key].fieldLineObj_f, rvs4)

    forward.DiffuseColor = [1, 0.0, 0.0]
    forward.LineWidth = 2.0


loc_spheres = {}
loc_disp = {}
for key in phi_location:
    print "Location: ", phi_location[key]
    loc_spheres[key] = pv.Sphere()
    loc_spheres[key].Radius = 0.01
    loc_spheres[key].Center = phi_location[key]
    loc_disp[key] = pv.Show(loc_spheres[key], rvs4)
    if np.any(np.isclose(key, calcTimes)):
        loc_disp[key].DiffuseColor = [0.0, 0.0, 1.0]
        loc_spheres[key].Radius = 0.011

# # create a new 'Poly Line Source'
# polyLineSource1 = pv.PolyLineSource()
#
# # Properties modified on polyLineSource1
# polyLineSource1.Closed = 1
#
# # show data in view
# polyLineSource1Display = pv.Show(polyLineSource1, rvs4)
# # trace defaults for the display properties.
# polyLineSource1Display.AmbientColor = [0.16993972686350806, 0.16993972686350806, 0.16993972686350806]
# polyLineSource1Display.ColorArrayName = [None, '']
# polyLineSource1Display.DiffuseColor = [1.0, 1.0, 0.0]
# polyLineSource1Display.BackfaceDiffuseColor = [0.3799038681620508, 0.3799038681620508, 0.3799038681620508]
# polyLineSource1Display.GlyphType = 'Arrow'
# polyLineSource1Display.CubeAxesColor = [0.16993972686350806, 0.16993972686350806, 0.16993972686350806]
# polyLineSource1Display.SetScaleArray = [None, '']
# polyLineSource1Display.ScaleTransferFunction = 'PiecewiseFunction'
# polyLineSource1Display.OpacityArray = [None, '']
# polyLineSource1Display.OpacityTransferFunction = 'PiecewiseFunction'


# points = []
# orderedKeys = sorted(phi_location.keys())

# for key in orderedKeys:
#     points.append(phi_location[key])

# pointsFlat = np.array(points).flatten()
# polyLineSource1.Points = pointsFlat

# tube1 = pv.Tube(Input=polyLineSource1)
# tube1.Radius = 0.005
# tube1Display = pv.Show(tube1, rvs4)
# tube1Display.DiffuseColor = [1.0, 1.0, 0.0]


pv.WriteImage("out/png/lstar_test.png", view=rvs4)
