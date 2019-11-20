from mpi4py import MPI
from ghostpy import algorithms as gpa
from ghostpy.Invariants import LShell as ls
from ghostpy import plotting as gpp
from ghostpy.data import VtkData as vtkd

class LstarWorker(object):

    def __init__(self, comm):
        assert isinstance(comm,MPI.Comm)
        self.rank = comm.Get_rank()
        self.comm = comm
        self.data = {}
        self.lstar = {}
        self.tData = None
        self.fData = None
        print("Worker Running on process {}".format(self.rank))

    def __call__(self, *args, **kwargs):
        status = self.__await_init_from_master__()
        status = self.__process_list__()
        status = self.__return_data__()
        print ("Process {} has completed".format(self.rank))

    def __await_init_from_master__(self):
        print ("Waiting for Init Commands")
        self.tData = self.comm.recv(source=0)
        self.fData = self.comm.recv(source=0)
        print("Initialization Message Received on rank {}".format(self.rank))

        for key in self.fData:
            self.data[key] = vtkd.VtkData(filename=self.fData[key][0], vector=self.fData[key][1])
            print("Initialized Data for File: {} on  process {}".format(self.fData[key][0], self.rank))

    def __process_list__(self):
        # process the lists of points and files
        # process points:

        if self.tData is  not None and len(self.data) > 0:
            for point in self.tData:
                self.lstar[tuple(point)] = ls.LShell(start_loc=point[:-1], data=self.data[point[3]]).l_star(res=10000)
                print ("[{}]\t-\tL*({})\t=\t{}".format(self.rank, point[:-1], self.lstar[tuple(point)]))

        return True

    def __request_instructions__(self):
        pass

    def __return_data__(self):
        self.comm.send(self.lstar, dest=0)
        print ("Data from Processor {} returned to master.".format(self.rank))
        pass


