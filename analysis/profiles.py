from __future__ import absolute_import

import numpy as np

import ghostpy.Invariants.LShell as ls
import ghostpy.algorithms.convert as algx
import ghostpy.algorithms.common as algc

import ghostpy.data.VtkData as vdt
import cProfile, pstats, StringIO
import ghostpy.data.GpData as gpd
import ghostpy.plotting.FieldLinePlot as flplt
import ghostpy.Invariants.FieldLine as fl
import random as rnd
import csv
from parse import *

import matplotlib.pylab as plt

import glob


def pointTest(pos=None, data=None):

    for loc in pos:
        lsl = ls.LShell(start_loc=loc, data=data, save_lines=True, error_tol=2e-6)
        lstar = lsl.l_star(res=10000)
        print ("LStar: {}".format(lstar))


def convergeTest(pos=None, data=None):
    for loc in pos:
        lsl = ls.LShell(start_loc=loc, data=data, save_lines=True, error_tol=2e-6)
        lstar = lsl.l_star(res=10000)
        print ("LStar: {}".format(lstar))
        lsl.converge_lstar(tol=1e-3)
        lstar = lsl.l_star(res=1000)
        print ("LStar(C): {}".format(lstar))


def convergeTest2(pos=None, data=None):
    for loc in pos:
        lsl = ls.LShell(start_loc=loc, data=data, save_lines=True, error_tol=2e-6)
        lstar = lsl.l_star(res=10000)
        print ("LStar: {}".format(lstar))
        lsl.converge_p2(depth=3)
        lstar = lsl.l_star(res=10000)
        print("LStar(C): {}".format(lstar))

def process_stats(prof=None):
    assert prof is not None
    ls_prof = {}
    timeParseFormat = "{:s}{:d} function calls ({:d} primitive calls) in {:f} seconds"
    timeParseFormat2 = "{:s}{:d}function calls in {:d} seconds"

    timeparse = parse(timeParseFormat, prof[0])

    if timeparse is None:
        timeparse = parse(timeParseFormat2, prof[0])

    ls_prof['time'] = timeparse[-1]

    keysline = prof[5].strip()
    keys = keysline.split()
    count = 0
    for line in prof[6:]:
        linedata = line.split()
        if len(linedata) == 0:
            continue
        linedatnum = [float(a) for a in linedata[:-1]]
        d = dict()
        for idx in range(len(linedatnum)):
            d[keys[idx]] = linedatnum[idx]
        ls_prof[linedata[-1]] = d
        count += 1
    return ls_prof


outpath = "/Volumes/8TB Seagate/PhD Data/profiles/"
data = vdt.VtkData(filename="../unit_tests/test_data/WHIQuad.vts", vector='B')

ps = [[2.43877196,  0.31925336,  0.74814534]]

profile = cProfile.Profile()
profile.enable()
pointTest(pos=ps, data=data)
profile.disable()
s1 = StringIO.StringIO()
ps1 = pstats.Stats(profile, stream=s1)
ps1 = ps1.strip_dirs()
ps1.print_stats('LShell.py:')

prof = s1.getvalue()
prof_lines1 = prof.split("\n")
lprofile = process_stats(prof_lines1)

fig = plt.figure()
ax = fig.add_subplot(111)






#
# s2 = StringIO.StringIO()
# ps2 = pstats.Stats(profile, stream=s2)
# ps2.print_stats('FieldLine.py:')
#
# prof2 = s2.getvalue()
# prof_lines = prof.split('\n')
# print(prof2)
#









