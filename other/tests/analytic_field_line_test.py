import fieldLine as fl
import numpy as np
from prototype import driftShell as ds

from ghostpy.prototype import inv_common as ih

fl1 = fl.fieldLine(start=[8., 0, 0])

x_f = [a[0] for a in fl1.fieldLinePoints_f]
y_f = [a[1] for a in fl1.fieldLinePoints_f]
z_f = [a[2] for a in fl1.fieldLinePoints_f]

x_b = [a[0] for a in fl1.fieldLinePoints_b]
y_b = [a[1] for a in fl1.fieldLinePoints_b]
z_b = [a[2] for a in fl1.fieldLinePoints_b]

print "B(K=0): {} nT".format(fl1.get_B_mirror_for_K(K=0))

divs = 1000
calcTimes = np.linspace(0, 24, num=4, endpoint=False)

ds1 = ds.driftShell(local_times=calcTimes, start_line=[2.5, 0.0, 0.0], K=0, mode='analytic_K')
ideal_RE = ih.mag(ds1.field_lines[0].fieldLinePoints_f[0])
print "Ideal L*: {}".format(ideal_RE)

test_grids = np.arange(start=100, stop=100, step=50, dtype=np.int)
L_star = ds1.calculate_L_star(RE=1.0, lat_divs=divs, lon_divs=divs)
diff = ideal_RE-L_star
error = np.abs(diff/ideal_RE * 100)

print "L*: {}".format(L_star)
print "diff: {}".format(diff)
print "% Error: {}".format(error)



# mlab.plot3d(x_f, y_f, z_f)
# mlab.plot3d(x_b, y_b, z_b)
#
# mlab.show()