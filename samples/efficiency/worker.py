from mpi4py import MPI
from ghostpy.Invariants import LShell as ls
import numpy as np
import time

class LStarConvergeWorker(object):

    def __init__(self, comm):
        assert isinstance(comm, MPI.Comm)
        self.rank = comm.Get_rank()
        self.data = None
        self.vector = None
        self.pointlist = None
        self.conv_tol = None
        self.loc = []
        self.lstar = []
        self.clstar = []
        self.clines = []
        self.kval = []
        self.bmirr = []
        self.time = []

        print ("Worker Process is Running on Rank {}".format(self.rank))

    def __call__(self, conv_tol=1e-6, data=None, vector=None, pointlist=None):
        self.conv_tol = conv_tol
        self.data = data
        self.vector = vector
        self.pointlist = pointlist
        self.__process_list__()
        print ("Points on Rank {} have finished.".format(self.rank))

    def __process_list__(self):
        if self.pointlist is not None:
            count = 0
            for point in self.pointlist:
                self.loc.append(point)
                dtimeS = time.time()
                lshell = ls.LShell(start_loc=point, data=self.data, error_tol=1e-6)
                self.lstar.append(lshell.l_star(res=1000))
                print ("[{} - {}] - L*({}) = {}".format(self.rank, count, point, self.lstar[-1]))

                if self.conv_tol > 0:
                    # lshell.converge_lstar(tol=1e-3)
                    lshell.converge_p2(depth=2)
                dtimeE = time.time()
                self.clstar.append(lshell.l_star(res=1000))
                self.clines.append(lshell.get_number_of_traces())
                self.bmirr.append(lshell.b)
                self.time.append(dtimeE-dtimeS)

                print ("[{} - {}] - Converged L*({}) = {} with {} lines".format(self.rank, count, point, self.clstar[-1], self.clines[-1]))

                self.kval.append(lshell.k)
                count += 1


    def get_data(self):
        return np.array(self.loc), np.array(self.lstar), np.array(self.clstar), np.array(self.clines), np.array(self.kval), np.array(self.bmirr), np.array(self.time)

