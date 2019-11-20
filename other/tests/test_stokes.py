import matplotlib.pyplot as plt
import numpy as np

from prototype import inv_common as ih

num_lats = 25

lats = np.linspace(start=0, stop=75, num=num_lats, endpoint=False)
lstar = []
lstarA = []
for lat in lats:
    list = ih.get_lat_path(lat=lat, divs=200)

    flux = ih.stokes_integration_r_theta_phi(path=list, val_fun=ih.get_dipole_value)
    lstar.append(ih.L_star_from_flux(flux))
    lstarA.append(ih.analytic_dipole_L(lat=lat, r=1.0))


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(lstar, lstarA, 'k.-', label="Calculated $L^{*}$ v. Analytic $L^{*}$")
ax.set_ylabel("Analytic $L^{*}$")
ax.set_xlabel("Calculated $L^{*}$")
ax.set_title("Stoke's Method v. Analytic $L^{*}$ on a Dipole")
fig.savefig("out/error_analysis/stokes_test.pdf")

plt.show()



