from mpi4py import MPI
from ghostpy.Invariants import LShell as ls
import numpy as np

class LStarConvergeWorker(object):

    def __init__(self, comm):
        assert isinstance(comm, MPI.Comm)
        self.rank = comm.Get_rank()
        self.data = None
        self.vector = None
        self.pointlist = None
        self.conv_tol = None
        self.k=0
        self.loc = []
        self.lstar = []
        self.clstar = []
        self.clines = []
        print ("Worker Process is Running on Rank {}".format(self.rank))

    def __call__(self, conv_tol=1e-6, data=None, vector=None, pointlist=None, k=0):
        self.conv_tol = conv_tol
        self.data = data
        self.vector = vector
        self.pointlist = pointlist
        self.__process_list__()
        self.k = k

    def __process_list__(self):
        if self.pointlist is not None:
            for point in self.pointlist:
                self.loc.append(point)
                lshell = ls.LShell(start_loc=point, k=self.k, data=self.data)
                self.lstar.append(lshell.l_star(res=10000))
                print ("[{}] - L*({}) = {}".format(self.rank, point, self.lstar[-1]))

                lshell.converge_lstar(tol=self.conv_tol)
                self.clstar.append(lshell.l_star(res=10000))
                self.clines.append(lshell.get_number_of_traces())
                print ("[{}] - Converged L*({}) = {} with {} lines".format(self.rank, point, self.clstar[-1], self.clines[-1]))

    def get_data(self):
        return np.array(self.loc), np.array(self.lstar), np.array(self.clstar), np.array(self.clines)

