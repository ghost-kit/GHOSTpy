from mpi4py import MPI
from ghostpy.Invariants import LShell as ls
from ghostpy.data import VtkData as vtkd
import numpy as np
import re

class LstarMaster(object):

    def __init__(self, comm):
        assert isinstance(comm, MPI.Comm)
        rank = comm.Get_rank()
        print ("Master Process Started on Rank: {}".format(rank))

        # TODO: will need to re-think the following dictionaries to handle same points at different times....
        # TODO: Currently it will delete duplicate points....
        self.fileTimeList = dict()
        self.pointTimeList = []
        self.calculatedLstar = dict()
        self.vector = None
        self.comm = comm
        self.points_per_proc = None
        self.points_remain = None
        self.num_workers = None
        self.worker_responsibility = dict()
        self.lstar = {}
        self.data = {}
        self.rank = rank
        self.in_lstar = {}
        self.out_file = None

    def __call__(self, *args, **kwargs):
        print("Master Process Running")
        self.numprocs = self.comm.Get_size()
        self.num_workers = int(self.numprocs)
        print ("Number of Workers: {}".format(self.num_workers))

        self.__part_out_points__()

        if self.num_workers > 1:
            self.__initialize_workers__()

        self.__process_local_points__()

        if self.num_workers > 1:
            self.__collect_data__()

        self.__write_data__()

        print("Processing Complete.")

    def add_file_time_pair(self, ftpair=None):
        assert isinstance(ftpair, tuple) and len(ftpair) == 3
        self.fileTimeList[ftpair[1]] = (ftpair[0], ftpair[2])

    def add_point_time_pair(self, ptpair=None):
        self.pointTimeList.append(ptpair)

    def get_number_of_points(self):
        return len(self.pointTimeList)

    def get_number_of_files(self):
        return len(self.fileTimeList)

    def __collect_data__(self):
        for proc in range(self.num_workers):
            if proc == 0:
                continue
            self.in_lstar[proc] = self.comm.recv(source=proc)
        print("Data Collected")

    def __write_data__(self):
        print("Writing Data To File")
        ofile = open("outdata.csv", "w+")
        lines = self.__format_data_for_output__(self.lstar)
        for l in self.in_lstar:
            self.__format_data_for_output__(self.in_lstar[l], lines=lines)

        ofile.writelines(lines)

    def __format_data_for_output__(self, input_dict, lines=None):
        if lines is None:
            lines = ["File Name,\t vector,\t SM Point (X),\t SM Point (Y),\t SM Point (Z),\t L* \n"]

        for dat in input_dict:
            filename = self.fileTimeList[np.array(dat)[-1]][0]
            vector = self.fileTimeList[np.array(dat)[-1]][1]
            ddat = np.array(dat)[:-1]
            ddat = str(ddat)
            ddat = ddat.replace("[", "")
            ddat = ddat.replace("]", "")
            ddat = re.sub("\s+",",\t ", ddat.strip())
            d = input_dict[dat]
            dline = str(filename) + ",\t" + str(vector) + ",\t" + str(ddat) + ",\t" + str(d) + '\n'
            lines.append(dline)

        return lines

    def __initialize_workers__(self):
        print ("Initializing Worker Processes")
        for proc in self.worker_responsibility:
            if proc == 0:
                continue
            self.comm.send(self.worker_responsibility[proc], dest=int(proc))
            self.comm.send(self.fileTimeList, dest=int(proc))

    def __process_local_points__(self):
        if len(self.worker_responsibility[self.rank]) > 0:
            print ("Master Process Working on Remaining items")

            for key in self.fileTimeList:
                self.data[key] = vtkd.VtkData(filename=self.fileTimeList[key][0], vector=self.fileTimeList[key][1])

            for point in self.worker_responsibility[self.rank]:
                self.lstar[tuple(point)] = ls.LShell(start_loc=point[:-1], data=self.data[point[3]]).l_star(res=10000)
                print ("[{}]\t-\tL*({})\t=\t{}".format(self.rank, point[:-1], self.lstar[tuple(point)]))

    def __part_out_points__(self):
        if self.num_workers != 0:
            self.points_per_proc = float(self.get_number_of_points()) / self.num_workers
            self.points_remain = self.get_number_of_points() % self.num_workers
            print ("Number of Points per processor: {}".format(self.points_per_proc))

        else:
            self.points_per_proc = 0
            self.points_remain = self.get_number_of_points()
        # set up the sorter

        for proc in range(self.num_workers):
            self.worker_responsibility[proc] = []

        # parting out points
        parts = self.get_number_of_points() # - self.points_remain

        for idx in range(parts):
            procNum = idx % self.num_workers
            self.worker_responsibility[procNum].append(self.pointTimeList[idx])

            # print("Process: {}, idx: {}".format(procNum, idx))
        # for idx in range(self.points_remain):
        #     self.remainder_pool.append(self.pointTimeList[parts+idx])

