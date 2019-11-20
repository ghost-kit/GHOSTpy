import numpy as np
import paraview.simple as pv

from prototype import driftShell as ds

pv._DisableFirstRenderCameraReset()

checkTimes = np.linspace(0, 24, num=64, endpoint=False)
#calcTimes = [0.0, 6.0, 12.0,  18.0]
calcTimes = np.linspace(0,24, num=4, endpoint=False)

#t96 = pv.OpenDataFile("../out/fields/t96_dp0_dst50_5.vts")
#t96 = pv.OpenDataFile("../out/fields/dipole_dp0.vts")
t96 = pv.OpenDataFile("../out/fields/t96_dp20_dst50_small_grid.vts")

t96_K0_ds2 = ds.driftShell(PV_dataObject=t96, local_times=checkTimes, b_mirror=150.00, K=100, mode='K')
t96_K0_ds = ds.driftShell(PV_dataObject=t96, local_times=calcTimes, b_mirror=150.00, K=100, mode='K')

RE_el = 1.0

phi_location = {}
phi_location2 = {}
print "Locations: ", phi_location

interp_points = np.linspace(0,24,64, endpoint=False)

for x in interp_points:
        phi_location[x] = t96_K0_ds.get_new_phi_for_localtime(x, RE=RE_el)
        phi_location2[x] = t96_K0_ds2.get_new_phi_for_localtime(x, RE=RE_el)

# now to figure out how to make the actual grid.


rvs4 = pv.CreateRenderView()
rvs4.InteractionMode = '2D'
rvs4.CameraPosition = [0, 0, 5]
rvs4.CameraFocalPoint = [0, 0, 0]
rvs4.CameraViewUp = [0, 1, 0]
rvs4.CameraParallelScale = 1.25
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
rvs4.ViewSize = [1920*5, 1280*5]
rvs4.Background = [1.0, 1.0, 1.0]

# create a sphere to represent Earth
earth = pv.Sphere()
earth.Radius = RE_el
earth.ThetaResolution = 128
earth.PhiResolution = 128

earth_disp4 = pv.Show(earth, rvs4)
earth_disp4.DiffuseColor = [0.51, 0.51, 0.51]

# display check lines
for key in t96_K0_ds2.field_lines:
    if t96_K0_ds2.field_lines[key] is None:
        continue
    forward = pv.Show(t96_K0_ds2.field_lines[key].fieldLineObj_f, rvs4)
    forward.DiffuseColor = [1, 0.0, 0.0]
    forward.LineWidth = 2.0


grid = t96_K0_ds.build_polar_cap_grid(RE=RE_el, divs=24, max_lat_divs=10)

locs = []
for key in grid:
    locs += grid[key]


loc_spheres = []
loc_disp = []
for loc in locs:
    loc_spheres.append(pv.Sphere())
    loc_spheres[-1].Radius = 0.005
    loc_spheres[-1].Center = loc
    loc_disp.append(pv.Show(loc_spheres[-1], rvs4))


# create a new 'Poly Line Source'
polyLineSource1 = pv.PolyLineSource()
polyLineSource2 = pv.PolyLineSource()

# Properties modified on polyLineSource1
polyLineSource1.Closed = 1
polyLineSource2.Closed = 1

# show data in view
polyLineSource1Display = pv.Show(polyLineSource1, rvs4)
polyLineSource2Display = pv.Show(polyLineSource2, rvs4)

pv.Hide(polyLineSource2, rvs4)
pv.Hide(polyLineSource1, rvs4)

# trace defaults for the display properties.
polyLineSource1Display.AmbientColor = [0.16993972686350806, 0.16993972686350806, 0.16993972686350806]
polyLineSource1Display.ColorArrayName = [None, '']
polyLineSource1Display.DiffuseColor = [1.0, 1.0, 0.0]
polyLineSource1Display.BackfaceDiffuseColor = [0.3799038681620508, 0.3799038681620508, 0.3799038681620508]

polyLineSource2Display.AmbientColor = [0.16993972686350806, 0.16993972686350806, 0.16993972686350806]
polyLineSource2Display.ColorArrayName = [None, '']
polyLineSource2Display.DiffuseColor = [1.0, 1.0, 0.0]
polyLineSource2Display.BackfaceDiffuseColor = [0.3799038681620508, 0.3799038681620508, 0.3799038681620508]


points = []
orderedKeys = sorted(phi_location.keys())

points2 = []
orderedKeys2 = sorted(phi_location2.keys())

for key in orderedKeys:
    points.append(phi_location[key])

for key in orderedKeys2:
    points2.append(phi_location2[key])


pointsFlat = np.array(points).flatten()
pointsFlat2 = np.array(points2).flatten()

polyLineSource1.Points = pointsFlat
polyLineSource2.Points = pointsFlat2

tube1 = pv.Tube(Input=polyLineSource1)
tube1.Radius = 0.005
tube1Display = pv.Show(tube1, rvs4)
tube1Display.DiffuseColor = [1.0, 1.0, 0.0]

tube2 = pv.Tube(Input=polyLineSource2)
tube2.Radius = 0.005
tube2Display = pv.Show(tube2, rvs4)
tube2Display.DiffuseColor = [0.0, 0.0, 1.0]

pv.WriteImage("out/png/polar_cap_grid_t96_dp20.png", view=rvs4)
