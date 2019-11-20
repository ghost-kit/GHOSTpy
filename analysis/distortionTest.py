from __future__ import absolute_import

import numpy as np

import ghostpy.Invariants.LShell as ls
import ghostpy.algorithms.convert as algx
import ghostpy.algorithms.common as algc
import ghostpy.algorithms.FieldTracers as ft
import ghostpy.algorithms.DipoleField as df
import ghostpy.data.VtkData as vdt
import ghostpy.data.LfmVtkData as lvdt
import ghostpy.data.DipoleData as dpd
import ghostpy.data.GpData as gpd
import ghostpy.plotting.FieldLinePlot as flplt
import ghostpy.Invariants.FieldLine as fl
import csv

import matplotlib.pylab as plt

datal = vdt.VtkData(filename="/Users/jjm390/src/ghostpy/ghostpy/unit_tests/test_data/WHIDouble.vts", vector="B")
ls1 = ls.LShell(data=datal, start_loc=[-13.3, 0.0, 0.0], save_lines=True)
ls1.converge_p2(depth=4)

plot = flplt.FieldLinePlot()
plot.plot_shell_field_lines(ls1)
plot.plot_drift_boundary(ls1)
plot.savePDF("/Volumes/8TB Seagate/PhD Data/profiles/distort.pdf")

plt.show()
