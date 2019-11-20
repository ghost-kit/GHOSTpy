# This is just a test to see if I can get 2 pi out of a unit hemisphere integration

import numpy as np

from prototype import inv_common as ih

# calcTimes = np.linspace(0,24, num=4, endpoint=False)
# t96 = pv.OpenDataFile("../out/fields/dipole_dp0_small_grid.vts")
# t96_K0_ds = ds.driftShell(PV_dataObject=t96, local_times=calcTimes, b_mirror=200.00, K=0, mode='K')
#
# RE_el = 1.0

solution = ih.analytic_sphere_cap_area(0.0, 1)
summation = 0
summation2 = 0

latitude = 0.0

lat_divs = 25
lon_divs = 25

print "Starting Check with ", lat_divs*lon_divs, "cells."

while not np.isclose(solution, summation, rtol=0.0, atol=0.005) and summation < solution * 1.1:
    summation = 0
    summation2 = 0
    lat_divs *= 2
    lon_divs *= 2
    hem_cart, hem_sph = ih.build_spherical_cap_grid(RE=1, grid_type="both", divs=lon_divs, base_lat=latitude)

    # make dpdt pieces
    elements = {}
    elements2 = {}
    for key in hem_sph:
        elements[key] = []
        elements2[key] = []

    hkeys = hem_sph.keys()
    num_lon = len(hkeys)

    # print "number of longitude keys: ", num_lon

    for index_lon in range(num_lon):
        for index_lat in range(len(hem_sph[hkeys[index_lon]])-1):

            lon = index_lon
            lat = index_lat

            if lon+1 >= num_lon:
                lower_lon = num_lon - 1
                upper_lon = 0
                # print "Readjusted longitude Key: ", lower_lon, ":", upper_lon

            else:
                lower_lon = lon
                upper_lon = lon + 1

            lower_lat = lat
            upper_lat = lat + 1

            # print "lower longitude index: ", lower_lon
            # print "upper longitude index: ", upper_lon
            # print "lower latitude index: ", lower_lat
            # print "upper latitude index: ", upper_lat

            r1 = hem_sph[hkeys[lower_lon]][lower_lat][0]
            lat1 = hem_sph[hkeys[lower_lon]][lower_lat][2]
            lat2 = hem_sph[hkeys[lower_lon]][upper_lat][2]
            lat3 = hem_sph[hkeys[upper_lon]][upper_lat][2]
            lat4 = hem_sph[hkeys[upper_lon]][upper_lat][2]

            lon1 = hem_sph[hkeys[lower_lon]][lower_lat][1]
            lon2 = hem_sph[hkeys[upper_lon]][lower_lat][1]
            lon3 = hem_sph[hkeys[upper_lon]][lower_lat][1]
            lon4 = hem_sph[hkeys[upper_lon]][lower_lat][1]

            coords1 = hem_cart[hkeys[lower_lon]][lower_lat]
            coords2 = hem_cart[hkeys[lower_lon]][upper_lat]
            coords3 = hem_cart[hkeys[upper_lon]][lower_lat]
            coords4 = hem_cart[hkeys[upper_lon]][upper_lat]

            lat1 = (np.pi/2)-lat1
            lat2 = (np.pi/2)-lat2

            if lon1 < 0:
                lon1 += 2*np.pi

            if lon2 < 0:
                lon2 += 2*np.pi

            # dT needs to be based on co-latitude, not latitude
            dT1 = lat1 - lat2
            dT2 = lat3 - lat4

            if lon2 < lon1:
                # print "Updating dP lon2"
                dP1 = (lon2 + (2 * np.pi)) - lon1
                dP2 = (lon4 + (2 * np.pi)) - lon3
            else:
                dP1 = lon2 - lon1
                dP2 = lon4 - lon3

            tri1 = 0.5 * dT1 * dP1
            tri2 = 0.5 * dT2 * dP2

            tot_area = tri1 + tri2

            elements[hkeys[index_lon]].append(r1 ** 2 * np.cos(lat1) * dT1 * dP1)
            # print "SP: ", elements[hkeys[index_lon]][-1]

            sphere_int_meth = (tot_area)
            # print "GC: ", sphere_int_meth
            elements2[hkeys[index_lon]].append(sphere_int_meth)

    for key in elements:
        summation += np.sum(elements[key])
        summation2 += np.sum(elements2[key])

    print "Cells: {} Ideal: {} Calc1: {} Diff: {}".format(int(lat_divs*lon_divs), solution, summation, solution - summation)
    # print "Cells: {} Ideal: {} Calc2: {} Diff: {}".format(int(lat_divs * lon_divs), solution, summation2, solution - summation2)

print "Final Lat_Divs: ", lat_divs
print "Final Lon_divs: ", lon_divs
print "Final Cell Count: ", lon_divs * lon_divs
print "Convergence Test Complete."

# phi_location = {}
# phi_location2 = {}
# print "Locations: ", phi_location
#
# interp_points = np.linspace(0,24,64, endpoint=False)
#
# for x in interp_points:
#     phi_location[x] = t96_K0_ds.get_new_phi_for_localtime(x, RE=RE_el)
#
# # now to figure out how to make the actual grid.
#
# rvs4 = pv.CreateRenderView()
# rvs4.InteractionMode = '2D'
# rvs4.CameraPosition = [0, 0, 5]
# rvs4.CameraFocalPoint = [0, 0, 0]
# rvs4.CameraViewUp = [0, 1, 0]
# rvs4.CameraParallelScale = 1.25
# rvs4.AxesGrid.Visibility = 1
# #rvs4.AxesGrid.AxesToLabel = 5
# rvs4.AxesGrid.GridColor = [0.8, 0.8, 0.8]
# rvs4.AxesGrid.ShowGrid = 1
# rvs4.AxesGrid.XLabelColor = [0.2, 0.2, 0.2]
# rvs4.AxesGrid.XTitleColor = [0.2, 0.2, 0.2]
# rvs4.AxesGrid.ZLabelColor = [0.2, 0.2, 0.2]
# rvs4.AxesGrid.ZTitleColor = [0.2, 0.2, 0.2]
# rvs4.AxesGrid.YLabelColor = [0.2, 0.2, 0.2]
# rvs4.AxesGrid.YTitleColor = [0.2, 0.2, 0.2]
# rvs4.OrientationAxesVisibility = 0
# rvs4.ViewSize = [1920*5, 1280*5]
# rvs4.Background = [1.0, 1.0, 1.0]
#
# # create a sphere to represent Earth
# earth = pv.Sphere()
# earth.Radius = RE_el
# earth.ThetaResolution = 128
# earth.PhiResolution = 128
#
# earth_disp4 = pv.Show(earth, rvs4)
# earth_disp4.DiffuseColor = [0.51, 0.51, 0.51]
#
#
# grid = hem_cart
#
# locs = []
# for key in grid:
#     locs += grid[key]
#
#
# loc_spheres = []
# loc_disp = []
# for loc in locs:
#     loc_spheres.append(pv.Sphere())
#     loc_spheres[-1].Radius = 0.005
#     loc_spheres[-1].Center = loc
#     loc_disp.append(pv.Show(loc_spheres[-1], rvs4))
#
# pv.WriteImage("out/png/hem_cart_grid.png", view=rvs4)
