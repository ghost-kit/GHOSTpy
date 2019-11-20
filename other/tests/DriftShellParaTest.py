import paraview.simple as pv
from mpi4py import MPI

from prototype import driftShellParalell as ds

# MPI stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# open a data set
# All ranks need the dataset, so lets open it on all ranks.
data = pv.OpenDataFile("../out/fields/t96_dp0_dst50_small_grid.vts")

#build the parallel drift shell
ds1 = ds.driftShellParalell(data, [0,6,12,18], K=0, b_mirror=115.0, mpi_comm=comm)

if rank is 0:
    print "Rank: ", rank, ": The only valid ds1 is on Rank 0, so lets deal with it here"
    print ds1.field_lines

else:
    print "Done with primary"



