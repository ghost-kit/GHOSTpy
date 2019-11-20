from __future__ import print_function

import sys
import numpy as np
import pickle
import ghostpy.data.VtkData as vtd
import ghostpy.Invariants.LShell as ls
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    print ("Output is being re-directed to rank-dependent files.  Please see these files for output.")

sys.stdout = open("stdout_rank{}.txt".format(rank), "w+", buffering=1)
sys.stderr = open("stderr_rank{}.txt".format(rank), "w+", buffering=1)

print("rank {}".format(rank))
print("rank {}".format(rank), file=sys.stderr)

pointMap = pickle.load(open("pointMap.pickle", "rb"))

numFiles = len(pointMap)

print ("Number of Files to Process: {}".format(numFiles))

workset = {}
ct = 0
for key in pointMap:
    if ct % size == rank:
        workset[key] = pointMap[key]
        print (key)
    ct += 1

print ("Rank: {} -- Files to Process: {}".format(rank, len(workset)))

RBSPA = {}
RBSPB = {}
ct = 1
for f in workset:
    print ("Rank: {} -- Processing {} of {} files".format(rank, ct, len(workset)))
    if rank == 0:
        print ("processing file: {}".format(f))
    data = vtd.VtkData(filename=f, vector="B")
    for pt in pointMap[f]:
        probe, tdiff, t, loc, alpha = pt
        loc = np.array(loc)
        lsp = ls.LShell(data=data, start_loc=loc, alpha=None)
        four_line = lsp.l_star(res=1000)
        print ("4-Point Lstar: {}".format(four_line))
        lsp.converge_lstar(tol=1e-3, int_res=1000)
        conv_line = lsp.l_star(res=1000)
        print ("Converged Lstar: {}".format(conv_line))
        k = lsp.k
        bm = lsp.b

        if probe == "A":
            print ("Saving to Probe A")
            RBSPA[t] = [four_line, conv_line, k, bm]
        else:
            print ("Saving to Probe B")
            RBSPB[t] = [four_line, conv_line, k, bm]
    ct += 1

RBSPA_dat = comm.gather(RBSPA, root=0)
RBSPB_dat = comm.gather(RBSPB, root=0)

if rank == 0:
    print ("Dumping to File")
    pickle.dump(RBSPA_dat, open("RBSPA_LFM_LKB.pickle", "wb"))
    pickle.dump(RBSPB_dat, open("RBSPB_LFM_LKB.pickle", "wb"))