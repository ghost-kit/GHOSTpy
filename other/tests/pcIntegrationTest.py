import numpy as np
import paraview.simple as pv
from prototype import driftShell as ds

from ghostpy.prototype import inv_common as ih

#pv._DisableFirstRenderCameraReset()

checkTimes = np.linspace(0, 24, num=64, endpoint=False)
calcTimes = np.linspace(0,24, num=4, endpoint=False)

t96 = pv.OpenDataFile("../out/fields/dipole_dp0_3_small_grid.vts")

t96_K0_ds = ds.driftShell(PV_dataObject=t96, local_times=calcTimes, b_mirror=115.00, K=0, mode='K')

RE_el = 1.0

L_10 = t96_K0_ds.calculate_L_star(RE=RE_el, lat_divs=10, lon_divs=10)
print "L* (10): ", L_10
raw_input("Press Enter to Continue...")

L_25 = t96_K0_ds.calculate_L_star(RE=RE_el, lat_divs=25, lon_divs=25)
print "L* (25): ", L_25
print "L* (10) - L* (25):", L_10 - L_25
raw_input("Press Enter to Continue...")

L_50 = t96_K0_ds.calculate_L_star(RE=RE_el, lat_divs=50, lon_divs=50)
print "L* (50): ", L_50
print "L* (25) - L* (50):", L_25-L_50
print "L* (10) - L* (50):", L_10-L_50
raw_input("Press Enter to Continue...")

L_100 = t96_K0_ds.calculate_L_star(RE=RE_el, lat_divs=100, lon_divs=100)
print "L* (100): ", L_100
print "L* (50) - L* (100):", L_50-L_100
print "L* (25) - L* (100):", L_25-L_100
print "L* (10) - L* (100):", L_10-L_100
raw_input("Press Enter to Continue...")

L_500 = t96_K0_ds.calculate_L_star(RE=RE_el, lat_divs=500, lon_divs=500)
print "L* (500): ", L_500
print "L* (100): ", L_100
print "L* (50): ", L_50
print "L* (25): ", L_25
print "L* (10): ", L_10
raw_input("Press Enter to Continue...")

#L_1000 = t96_K0_ds.calculate_L_star(RE=RE_el, lat_divs=1000, lon_divs=1000)
#print "L* (1000): ", L_1000
# raw_input("Press Enter to Continue...")

rvs4 = pv.CreateRenderView()
rvs4.InteractionMode = '2D'
rvs4.CameraPosition = [0, -190, 0]
rvs4.CameraFocalPoint = [0, 0, 0]
rvs4.CameraViewUp = [0, 0, 1]
rvs4.CameraParallelScale = 50
rvs4.AxesGrid.Visibility = 1
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

L_Text = pv.Text()
L_Text.Text = "$L^*=$ " + str(L_100)
Ldisp = pv.Show(L_Text, rvs4)
Ldisp.WindowLocation = 'LowerCenter'


for key in t96_K0_ds.field_lines:
    if t96_K0_ds.field_lines[key] is None:
        continue
    forward = pv.Show(t96_K0_ds.field_lines[key].fieldLineObj_f, rvs4)
    backward = pv.Show(t96_K0_ds.field_lines[key].fieldLineObj_b, rvs4)

    forward.DiffuseColor = [1, 0.0, 0.0]
    forward.LineWidth = 2.0
    backward.DiffuseColor = [1, 0.0, 0.0]
    backward.LineWidth = 2.0


eqc = ih.mag(t96_K0_ds.field_lines[0.0].fieldLinePoints_f[0])


print "Equatorial Crossing of Drift Shell: ", eqc
print "L* (500): ", L_500
print "L* (100): ", L_100
print "L* (50): ", L_50
print "L* (25): ", L_25
print "L* (10): ", L_10
print
print "Error: "
print "-------"
print "div=10: ", (eqc - L_10)/eqc * 100, "%"
print "div=25: ", (eqc - L_25)/eqc * 100, "%"
print "div=50: ", (eqc - L_50)/eqc * 100, "%"
print "div=100: ", (eqc - L_100)/eqc * 100, "%"
print "div=500: ", (eqc - L_500)/eqc * 100, "%"

pv.RenderAllViews()
pv.WriteImage("out/png/lstar_test.png", view=rvs4)
