import numpy as np
import sys
from mpi4py import MPI

import driftShell as ds


class driftShellParalell(ds.driftShell):

    def __init__(self, PV_dataObject, local_times=None, K=None, b_mirror=None, mode='K', mpi_comm=None):
        """
        This is the parallel version of the driftShell routines
        :param PV_dataObject: ParaView data object that contians the grid and data
        :param local_times: the seed local times for drift shell calculation
        :param K: value of 'K' used to calculate the drift shell
        :param b_mirror: B value used to cacluculate the drift shell
        :param mode: Mode to run the system in.  Please use 'K' mode, as 'location' has not been finished yet.
        """

        if mpi_comm is None:
            self.comm = MPI.COMM_WORLD
        else:
            self.comm = mpi_comm

        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        self.data = PV_dataObject
        self.is_valid = True

        if mode is 'K':
            # collective location
            self.field_lines = dict()
            self.start_RE = 5.0
            self.b_mirror = b_mirror
            self.K = K

            if self.rank == 0:
                print "Running in 'K' based mode"

            if self.rank == 0:
                workerCount = self.size-1
                if len(local_times) > workerCount:
                    jobs_per_worker = len(local_times) % workerCount
                else:
                    jobs_per_worker = 0

                print "Controller", self.rank
                print "Number of workers: ", workerCount
                print "Jobs per worker: ", jobs_per_worker

                if jobs_per_worker is 0:
                    print "Sending single job to each processor until jobs are done, then sending None"
                    lt = np.full((self.size, len(local_times)), -1, dtype=np.int)

                    for dest in range(1, workerCount):
                        # check to see if we are in range
                        if dest-1 < len(local_times):
                            lt[dest,0] = local_times[dest-1]
                    print "LT: ", lt

            else:
                lt = np.empty((self.size, len(local_times)), dtype=np.int)

            #collective broadcast
            self.comm.Bcast(lt, root=0)

            if self.rank is not 0:
                print "Rank: ", self.rank, ": Working Set: ", lt[self.rank]
                for x in lt[self.rank]:
                    if x >= 0:
                        print "Processing Local Time: ", x
                        self.add_field_line_from_K(local_time=x)

                    else:
                        print "No more work. "
                        break

            print "Gathering on Rank: ", self.rank
            self.field_lines = self.comm.gather(self.field_lines, root=0)


        else:
            print "Mode ", mode, " is not supported. Exiting"
            sys.exit()



