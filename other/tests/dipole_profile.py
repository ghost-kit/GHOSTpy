import collections as col
import matplotlib.pylab as plt
import numpy as np

from prototype import inv_common as ic

lats = np.linspace(start=0, stop=75, num=2500, endpoint=False)

points = col.OrderedDict()

for lat in lats:
    points[lat] = ic.analytic_dipole_L(lat=lat, r=1)

Y,X = ic.dict_to_x_y(points)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("Analytic $L^{*}$")
ax.set_ylabel("Intersection $\lambda$")
ax.set_title("Analytic $L^{*}$ vs. Intersection $\lambda$ Profile\non a Dipole")
ax.grid()
ax.plot(X, Y)

fig.savefig("out/profiles/dipole_Lstar_lambda_profile.pdf")

