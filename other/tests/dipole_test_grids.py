import numpy as np

from pyevtk.hl import gridToVTK

from prototype import bfieldgrid

re_res = np.linspace(start=1, stop=19, num=10, dtype=int)

for res in re_res:
    print "Creating Grid of Resolution: {}RE".format(res)
    E = bfieldgrid.bFieldGrid(extents=[-15, 15, -15, 15, -15, 15], resolution=(31 * res, 31 * res, 31 * res), exclusion_distance=0.2)
    gridToVTK("../../out/fields/dipole_dp0_{}_grid".format(res), E.x, E.y, E.z, pointData={"B": E.dipole_field(ps=0.0)})

