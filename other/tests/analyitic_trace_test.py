# testing field line trace

# import paraview.simple as pv
import matplotlib.pyplot as plt
import numpy as np

import mayavi.mlab as mav

from prototype import inv_common as ih

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.grid(True)
ax.set_xlim((-8, 8))
ax.set_ylim((-8, 8))
ax.set_zlim((-8, 8))

circ = plt.Circle((0, 0), 1, color='k', fill=False)
ax.add_artist(circ)
ax.auto_scale_xyz([-8, 8], [-8, 8], [-8, 8])

s_u = np.linspace(0, 2 * np.pi, 100)
s_v = np.linspace(0, np.pi, 100)

s_x = np.outer(np.cos(s_u), np.sin(s_v))
s_y = np.outer(np.sin(s_u), np.sin(s_v))
s_z = np.outer(np.ones(np.size(s_u)), np.cos(s_v))

earth = mav.mesh(s_x, s_y, s_z)

lts = np.linspace(0, 24, 24, endpoint=False)
loc = []
for lt in lts:
    loc.append(ih.get_location_from_localtime(lt)*6.0)

flines = {}
blines = {}
for line in loc:
    # for line in [6.0, -6.0]:
    print "Calculating Line: {}".format(line)
    x0 = line

    x = x0
    dpx = np.array(ih.get_dipole_value([x])[0])
    dt = 1e-7
    ib = 0.90
    ms = 1e6

    path_f, value_f = ih.trace_rk45(inner_boundary=ib, h=dt, x0=x0, max_steps=int(ms), direction='f')
    path_b, value_b = ih.trace_rk45(inner_boundary=ib, h=dt, x0=x0, max_steps=int(ms), direction='b')

    X_f = [a[0] for a in path_f]
    Y_f = [a[1] for a in path_f]
    Z_f = [a[2] for a in path_f]

    X_b = [a[0] for a in path_b]
    Y_b = [a[1] for a in path_b]
    Z_b = [a[2] for a in path_b]

    flines[tuple(line)] = [X_f, Y_f, Z_f]
    blines[tuple(line)] = [X_b, Y_b, Z_b]

    ax.plot(X_f, Y_f, Z_f, color='k')
    ax.plot(X_b, Y_b, Z_b, color='b')
    # ax.plot(X_b, Y_b, color='r')
    lf = mav.plot3d(X_f, Y_f, Z_f)
    lb = mav.plot3d(X_b, Y_b, Z_b)


# plt.show()
# mav.view(42, 73, 104, [0,  0,  0])

mav.show()
