import matplotlib.pyplot as plt
import numpy as np

from prototype import inv_common as ih

lats = np.arange(start=0, stop=90, step=0.5)
divs = 50
data = {}
data2 = {}
error = {}
error2 = {}

for lat in lats:
    grid_c, grid_s = ih.build_spherical_cap_grid(RE=1, divs=divs, base_lat=lat, grid_type='both')
    keys = grid_s.keys()
    r, theta, phi = grid_s[keys[0]][0]
    x, y, z = grid_c[keys[0]][0]
    dat2 = ih.get_lat((x, y, z))
    data[lat] = ih.rad_to_deg(theta)
    data2[lat] = dat2
    error[lat] = np.abs(lat - data[lat])
    error2[lat] = np.abs(lat - dat2)
    print (x, y, z)

    print ("Lat = {}, Grid Lat = {}".format(lat, dat2))

X, Y = ih.dict_to_x_y(data)

Xe, Ye = ih.dict_to_x_y(error)
Xe2, Ye2 = ih.dict_to_x_y(error2)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
line1, = ax1.plot(X,Y, 'r-', label="correlation of latitude")

ax2 = ax1.twinx()
line2, = ax2.plot(Xe,Ye, 'k-', label="Absolute Error (sphere)")
line3, = ax2.plot(Xe2, Ye2, 'b--', label="Absolute Error (cart)")

ax1.set_ylabel("Returned $\lambda$")
ax1.set_xlabel("Requested $\lambda$")
ax2.set_ylabel("Absolute Error")
ax1.set_title("Grid Generation Errors by Starting $\lambda$")

plt.legend(handles=[line1, line2, line3], loc=2)
plt.show()
