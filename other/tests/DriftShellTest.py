# Test drift shells routines

import paraview.simple as pv

from prototype import driftShell as ds

pv._DisableFirstRenderCameraReset()

# Load Test Data
t96 = pv.OpenDataFile("../out/fields/t96_128_dst50_5.vts")
dipole = pv.OpenDataFile("../out/fields/dipole_5.vts")

t96_K0_ds = ds.driftShell(PV_dataObject=t96, local_times=[12, 0], b_mirror=115.0, K=0, mode='K')

b_k1000 = t96_K0_ds.field_lines[12].get_B_mirror_for_K(K=1000)
t96_K1000_ds = ds.driftShell(PV_dataObject=t96, local_times=[12, 0], b_mirror=b_k1000, K=1000, mode='K')

dipole_K0_ds = ds.driftShell(PV_dataObject=dipole, local_times=[12, 0], b_mirror=115.0, K=0, mode='K')

b_k1000_d = dipole_K0_ds.field_lines[12].get_B_mirror_for_K(K=1000)
dipole_K1000_ds = ds.driftShell(PV_dataObject=dipole, local_times=[12, 0], b_mirror=b_k1000_d, K=1000, mode='K')

t96_K0_ds_m = ds.driftShell(PV_dataObject=t96, local_times=[0, 12], b_mirror=115.0, K=0, mode='K')

b_k1000_m = t96_K0_ds_m.field_lines[0].get_B_mirror_for_K(K=1000)
t96_K1000_midnight_ds = ds.driftShell(PV_dataObject=t96, local_times=[0, 12], b_mirror=b_k1000_m, K=1000, mode='K')



# print "Drift Shell [B_mirror(K)] (Dipole K0)"
# for localTime in dipole_K0_ds.field_lines:
#     print "B(K0): ", dipole_K0_ds.field_lines[localTime].get_B_mirror_for_K(0)
#     if dipole_K0_ds.field_lines[localTime] is None:
#         print "LT: ", localTime, ": NO LINE FOUND"
#     else:
#         print "LT: ", localTime, ": ", dipole_K0_ds.field_lines[localTime].startLoc
#
# print "Drift Shell [B_mirror(K)] (Dipole K1000)"
# for localTime in dipole_K1000_ds.field_lines:
#     print "B(K1000): ", dipole_K1000_ds.field_lines[localTime].get_B_mirror_for_K(1000)
#     if dipole_K1000_ds.field_lines[localTime] is None:
#         print "LT: ", localTime, ": NO LINE FOUND"
#     else:
#         print "LT: ", localTime, ": ", dipole_K1000_ds.field_lines[localTime].startLoc
#
# print "Drift Shell [B_mirror(K)] (T96 K0)"
# for localTime in t96_K0_ds.field_lines:
#     print "B(K0): ", t96_K0_ds.field_lines[localTime].get_B_mirror_for_K(0)
#     if t96_K0_ds.field_lines[localTime] is None:
#         print "LT: ", localTime, ": NO LINE FOUND"
#     else:
#         print "LT: ", localTime, ": ", t96_K0_ds.field_lines[localTime].startLoc
#
# print "Drift Shell [B_mirror(K)] (T96 K1000)"
# for localTime in t96_K1000_ds.field_lines:
#     print "B(K1000): ", t96_K1000_ds.field_lines[localTime].get_B_mirror_for_K(1000)
#     if t96_K1000_ds.field_lines[localTime] is None:
#         print "LT: ", localTime, ": NO LINE FOUND"
#     else:
#         print "LT: ", localTime, ": ", t96_K1000_ds.field_lines[localTime].startLoc
#
#


# Create a paraview render view so we can see visual progress.
# pv._DisableFirstRenderCameraReset()
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
rvs.AxesGrid.YLabelColor = [0.2, 0.2, 0.2]
rvs.AxesGrid.YTitleColor = [0.2, 0.2, 0.2]
rvs.OrientationAxesVisibility = 0
rvs.ViewSize = [1920, 1280]
rvs.Background = [1.0, 1.0, 1.0]


rvs2 = pv.CreateRenderView()
rvs2.InteractionMode = '2D'
rvs2.CameraPosition = [0, -5, 0]
rvs2.CameraFocalPoint = [0, 0, 0]
rvs2.CameraViewUp = [0, 0, 1]
rvs2.CameraParallelScale = 11
rvs2.AxesGrid.Visibility = 1
rvs2.AxesGrid.AxesToLabel = 5
rvs2.AxesGrid.GridColor = [0.8, 0.8, 0.8]
rvs2.AxesGrid.ShowGrid = 1
rvs2.AxesGrid.XLabelColor = [0.2, 0.2, 0.2]
rvs2.AxesGrid.XTitleColor = [0.2, 0.2, 0.2]
rvs2.AxesGrid.ZLabelColor = [0.2, 0.2, 0.2]
rvs2.AxesGrid.ZTitleColor = [0.2, 0.2, 0.2]
rvs2.AxesGrid.YLabelColor = [0.2, 0.2, 0.2]
rvs2.AxesGrid.YTitleColor = [0.2, 0.2, 0.2]
rvs2.OrientationAxesVisibility = 0
rvs2.ViewSize = [1920, 1280]
rvs2.Background = [1.0, 1.0, 1.0]

rvs3 = pv.CreateRenderView()
rvs3.InteractionMode = '2D'
rvs3.CameraPosition = [0, -5, 0]
rvs3.CameraFocalPoint = [0, 0, 0]
rvs3.CameraViewUp = [0, 0, 1]
rvs3.CameraParallelScale = 11
rvs3.AxesGrid.Visibility = 1
rvs3.AxesGrid.AxesToLabel = 5
rvs3.AxesGrid.GridColor = [0.8, 0.8, 0.8]
rvs3.AxesGrid.ShowGrid = 1
rvs3.AxesGrid.XLabelColor = [0.2, 0.2, 0.2]
rvs3.AxesGrid.XTitleColor = [0.2, 0.2, 0.2]
rvs3.AxesGrid.ZLabelColor = [0.2, 0.2, 0.2]
rvs3.AxesGrid.ZTitleColor = [0.2, 0.2, 0.2]
rvs3.AxesGrid.YLabelColor = [0.2, 0.2, 0.2]
rvs3.AxesGrid.YTitleColor = [0.2, 0.2, 0.2]
rvs3.OrientationAxesVisibility = 0
rvs3.ViewSize = [1920, 1280]
rvs3.Background = [1.0, 1.0, 1.0]

rvs4 = pv.CreateRenderView()
rvs4.InteractionMode = '2D'
rvs4.CameraPosition = [0, 0, 5]
rvs4.CameraFocalPoint = [0, 0, 0]
rvs4.CameraViewUp = [0, 1, 0]
rvs4.CameraParallelScale = 11
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
earth_disp = pv.Show(earth, rvs)
earth_disp2 = pv.Show(earth, rvs2)
earth_disp3 = pv.Show(earth, rvs3)
earth_disp4 = pv.Show(earth, rvs4)
earth_disp4.DiffuseColor = [0.81, 0.81, 0.81]
earth_disp3.DiffuseColor = [0.333, 0.333, 0.333]
earth_disp2.DiffuseColor = [0.333, 0.333, 0.333]
earth_disp.DiffuseColor = [0.333, 0.333, 0.333]



Bmag_calc = pv.Calculator(Input=t96)
Bmag_calc.ResultArrayName = 'Bmag'
Bmag_calc.Function ='mag(B)'
Bmag_slice = pv.Slice(Input=Bmag_calc)
Bmag_slice.SliceType = 'Plane'
Bmag_slice.SliceOffsetValues = [0.0]
Bmag_slice.SliceType.Origin = [-20.0, 0.0, 0.0]
Bmag_slice.SliceType.Normal = [0.0, 1.0, 0.0]

Bmag_calc2 = pv.Calculator(Input=dipole)
Bmag_calc2.ResultArrayName = 'Bmag'
Bmag_calc2.Function ='mag(B)'
Bmag_slice2 = pv.Slice(Input=Bmag_calc2)
Bmag_slice2.SliceType = 'Plane'
Bmag_slice2.SliceOffsetValues = [0.0]
Bmag_slice2.SliceType.Origin = [-20.0, 0.0, 0.0]
Bmag_slice2.SliceType.Normal = [0.0, 1.0, 0.0]


contour1 = pv.Contour(Input=Bmag_slice)
contour1.ContourBy = ['POINTS','Bmag']
contour1.Isosurfaces=[115, b_k1000]
cont1disp = pv.Show(contour1, rvs)
cont1disp.DiffuseColor = [0.0, 0.0, 0.0]

contour2 = pv.Contour(Input=Bmag_slice2)
contour2.ContourBy = ['POINTS','Bmag']
contour2.Isosurfaces=[115, b_k1000_d]
cont2disp = pv.Show(contour2, rvs2)
cont2disp.DiffuseColor = [0.0, 0.0, 0.0]

contour3 = pv.Contour(Input=Bmag_slice)
contour3.ContourBy = ['POINTS','Bmag']
contour3.Isosurfaces=[115, b_k1000_m]
cont1disp = pv.Show(contour3, rvs3)
cont1disp.DiffuseColor = [0.0, 0.0, 0.0]


b_loc_dipole_K0 = {}
b_loc_dipole_K1000 = {}
b_loc_t96_K0 = {}
b_loc_t96_K1000 = {}
b_loc_t96_K0_m = {}
b_loc_t96_K1000_m = {}

t96_k0_spheres = {}
t96_k0_spheres2 = {}
t96_k0_spheres3 = {}
t96_k0_spheres4 = {}

t96_k1000_spheres = {}
t96_k1000_spheres2 = {}
t96_k1000_spheres3 = {}
t96_k1000_spheres4 = {}

dipole_k0_spheres = {}
dipole_k0_spheres2 = {}


dipole_K1000_spheres = {}
dipole_K1000_spheres2 = {}


for lt in dipole_K0_ds.field_lines:
    print "Local Time: ", lt
    b_loc_dipole_K0[lt] = dipole_K0_ds.field_lines[lt].get_location_for_B(115)
    print "Locations: ", lt, ": ", b_loc_dipole_K0


for lt in dipole_K1000_ds.field_lines:
    b_loc_dipole_K1000[lt] = dipole_K1000_ds.field_lines[lt].get_location_for_B(b_k1000_d)

for lt in t96_K0_ds.field_lines:
    b_loc_t96_K0[lt] = t96_K0_ds.field_lines[lt].get_location_for_B(115)

for lt in t96_K1000_ds.field_lines:
    b_loc_t96_K1000[lt] = t96_K1000_ds.field_lines[lt].get_location_for_B(b_k1000)

for lt in t96_K0_ds_m.field_lines:
    b_loc_t96_K0_m[lt] = t96_K0_ds_m.field_lines[lt].get_location_for_B(115)

for lt in t96_K1000_midnight_ds.field_lines:
    b_loc_t96_K1000_m[lt] = t96_K1000_midnight_ds.field_lines[lt].get_location_for_B(b_k1000_m)

for key in t96_K0_ds.field_lines:
    if t96_K0_ds.field_lines[key] is None:
        continue
    forward = pv.Show(t96_K0_ds.field_lines[key].fieldLineObj_f, rvs)

    backward = pv.Show(t96_K0_ds.field_lines[key].fieldLineObj_b, rvs)

    forward.DiffuseColor = [0.6666666666666666, 0.0, 0.0]
    backward.DiffuseColor = [0., 0., 0.66]

    t96_k0_spheres[key] = pv.Sphere()
    t96_k0_spheres[key].Center = b_loc_t96_K0[key][0]
    t96_k0_spheres[key].Radius = 0.1
    t96_k0_spheres[key].ThetaResolution = 4
    t96_k0_spheres[key].PhiResolution = 4
    disp = pv.Show(t96_k0_spheres[key], rvs)
    disp.DiffuseColor = [.5, .5, .5]

    t96_k0_spheres2[key] = pv.Sphere()
    t96_k0_spheres2[key].Center = b_loc_t96_K0[key][1]
    t96_k0_spheres2[key].Radius = 0.1
    t96_k0_spheres2[key].ThetaResolution = 4
    t96_k0_spheres2[key].PhiResolution = 4
    disp = pv.Show(t96_k0_spheres2[key], rvs)

    disp.DiffuseColor = [.5, .5, .5]

for key in t96_K1000_ds.field_lines:
    if t96_K1000_ds.field_lines[key] is None:
        continue
    forward2 = pv.Show(t96_K1000_ds.field_lines[key].fieldLineObj_f, rvs)
    backward2 = pv.Show(t96_K1000_ds.field_lines[key].fieldLineObj_b, rvs)
    forward2.DiffuseColor = [0.6666666666666666, 0.0, 0.0]
    backward2.DiffuseColor = [0., 0., 0.66]

    t96_k1000_spheres[key] = pv.Sphere()
    t96_k1000_spheres[key].Center = b_loc_t96_K1000[key][0]
    t96_k1000_spheres[key].Radius = 0.1
    t96_k1000_spheres[key].ThetaResolution = 8
    t96_k1000_spheres[key].PhiResolution = 8
    disp = pv.Show(t96_k1000_spheres[key], rvs)
    disp.DiffuseColor = [.75, .75, .75]

    t96_k1000_spheres2[key] = pv.Sphere()
    t96_k1000_spheres2[key].Center = b_loc_t96_K1000[key][1]
    t96_k1000_spheres2[key].Radius = 0.1
    t96_k1000_spheres2[key].ThetaResolution = 8
    t96_k1000_spheres2[key].PhiResolution = 8
    disp = pv.Show(t96_k1000_spheres2[key], rvs)
    disp.DiffuseColor = [.75, .75, .75]

for key in t96_K0_ds_m.field_lines:
    if t96_K0_ds_m.field_lines[key] is None:
        continue
    forward = pv.Show(t96_K0_ds_m.field_lines[key].fieldLineObj_f, rvs3)
    backward = pv.Show(t96_K0_ds_m.field_lines[key].fieldLineObj_b, rvs3)
    forward.DiffuseColor = [0.6666666666666666, 0.0, 0.0]
    backward.DiffuseColor = [0., 0., 0.66]

    t96_k0_spheres3[key] = pv.Sphere()
    t96_k0_spheres3[key].Center = b_loc_t96_K0_m[key][0]
    t96_k0_spheres3[key].Radius = 0.1
    t96_k0_spheres3[key].ThetaResolution = 4
    t96_k0_spheres3[key].PhiResolution = 4
    disp = pv.Show(t96_k0_spheres3[key], rvs3)
    disp.DiffuseColor = [.5, .5, .5]

    t96_k0_spheres4[key] = pv.Sphere()
    t96_k0_spheres4[key].Center = b_loc_t96_K0_m[key][1]
    t96_k0_spheres4[key].Radius = 0.1
    t96_k0_spheres4[key].ThetaResolution = 4
    t96_k0_spheres4[key].PhiResolution = 4
    disp = pv.Show(t96_k0_spheres4[key], rvs3)

    disp.DiffuseColor = [.5, .5, .5]

for key in t96_K1000_midnight_ds.field_lines:
    if t96_K1000_midnight_ds.field_lines[key] is None:
        continue
    forward2 = pv.Show(t96_K1000_midnight_ds.field_lines[key].fieldLineObj_f, rvs3)
    backward2 = pv.Show(t96_K1000_midnight_ds.field_lines[key].fieldLineObj_b, rvs3)
    forward2.DiffuseColor = [0.6666666666666666, 0.0, 0.0]
    backward2.DiffuseColor = [0., 0., 0.66]

    t96_k1000_spheres3[key] = pv.Sphere()
    t96_k1000_spheres3[key].Center = b_loc_t96_K1000_m[key][0]
    t96_k1000_spheres3[key].Radius = 0.1
    t96_k1000_spheres3[key].ThetaResolution = 8
    t96_k1000_spheres3[key].PhiResolution = 8
    disp = pv.Show(t96_k1000_spheres3[key], rvs3)
    disp.DiffuseColor = [.75, .75, .75]

    t96_k1000_spheres4[key] = pv.Sphere()
    t96_k1000_spheres4[key].Center = b_loc_t96_K1000_m[key][1]
    t96_k1000_spheres4[key].Radius = 0.1
    t96_k1000_spheres4[key].ThetaResolution = 8
    t96_k1000_spheres4[key].PhiResolution = 8
    disp = pv.Show(t96_k1000_spheres4[key], rvs3)
    disp.DiffuseColor = [.75, .75, .75]


for key in dipole_K0_ds.field_lines:
    if dipole_K0_ds.field_lines[key] is None:
        continue
    forward = pv.Show(dipole_K0_ds.field_lines[key].fieldLineObj_f, rvs2)
    backward = pv.Show(dipole_K0_ds.field_lines[key].fieldLineObj_b, rvs2)
    forward.DiffuseColor = [0.6666666666666666, 0.0, 0.0]
    backward.DiffuseColor = [0., 0., 0.66]

    dipole_k0_spheres[key] = pv.Sphere()
    dipole_k0_spheres[key].Center = b_loc_dipole_K0[key][0]
    dipole_k0_spheres[key].Radius = 0.1
    dipole_k0_spheres[key].ThetaResolution = 4
    dipole_k0_spheres[key].PhiResolution = 4
    disp = pv.Show(dipole_k0_spheres[key], rvs2)
    disp.DiffuseColor = [.5, .5, .5]

    dipole_k0_spheres2[key] = pv.Sphere()
    dipole_k0_spheres2[key].Center = b_loc_dipole_K0[key][1]
    dipole_k0_spheres2[key].Radius = 0.1
    dipole_k0_spheres2[key].ThetaResolution = 4
    dipole_k0_spheres2[key].PhiResolution = 4
    disp = pv.Show(dipole_k0_spheres2[key], rvs2)
    disp.DiffuseColor = [.5, .5, .5]


for key in dipole_K1000_ds.field_lines:
    if dipole_K1000_ds.field_lines[key] is None:
        continue
    forward2 = pv.Show(dipole_K1000_ds.field_lines[key].fieldLineObj_f, rvs2)
    backward2 = pv.Show(dipole_K1000_ds.field_lines[key].fieldLineObj_b, rvs2)
    forward2.DiffuseColor = [0.6666666666666666, 0.0, 0.0]
    backward2.DiffuseColor = [0., 0., 0.66]

    dipole_K1000_spheres[key] = pv.Sphere()
    dipole_K1000_spheres[key].Center = b_loc_dipole_K1000[key][0]
    dipole_K1000_spheres[key].Radius = 0.1
    dipole_K1000_spheres[key].ThetaResolution = 8
    dipole_K1000_spheres[key].PhiResolution = 8
    disp = pv.Show(dipole_K1000_spheres[key], rvs2)
    disp.DiffuseColor = [.75, .75, .75]

    dipole_K1000_spheres2[key] = pv.Sphere()
    dipole_K1000_spheres2[key].Center = b_loc_dipole_K1000[key][1]
    dipole_K1000_spheres2[key].Radius = 0.1
    dipole_K1000_spheres2[key].ThetaResolution = 8
    dipole_K1000_spheres2[key].PhiResolution = 8
    disp = pv.Show(dipole_K1000_spheres2[key], rvs2)
    disp.DiffuseColor = [.75, .75, .75]

pv.RenderAllViews()

pv.WriteImage("drift_shell_t96.png", view=rvs)
pv.WriteImage("drift_shell_dipole.png", view=rvs2)
pv.WriteImage("drift_shell_t96_m.png", view=rvs3)



