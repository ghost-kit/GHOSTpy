import numpy as np

from prototype import inv_common as ih

Ls = np.linspace(start=2, stop=8, num=4800)

for L in Ls:
    print  "Crossing for L={}:\t{}".format(L,ih.analytic_dipole_line_lat_crossing(L=L, r=1))

