# File: parallelLstar.py
# Author: Joshua J. Murphy
# Project: Dissertation analysis
# Purpose: Processes LStar points in parallel utilizing MPI
#
# =============================================================

from mpi4py import MPI
import sys

# setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def main():
    if rank == 0:
        print ("Output is being re-directed to rank-dependent files.  Please see these files for output.")
    sys.stdout = open("stdout_rank{}.txt".format(rank), "w+", buffering=1)
    sys.stderr = open("stderr_rank{}.txt".format(rank), "w+", buffering=1)
    if rank == 0:
        ###################
        #   MASTER RANK   #
        ###################

        # master dependent modules
        import lsArgs
        import master
        from ghostpy.data import VtkData as gpd
        from ghostpy.algorithms import common as algc
        import numpy as np

        # get the arguments from the command line
        usage = "System Usage: {} -f <filename> -v <vector> -p <pointsFile>, -r <pointsRadius>".format(sys.argv[0])
        valid_args = {'f': "filename", 'v': "vector", 'r': "radius", 'p': "pointsFile"}
        arg_dict = lsArgs.lsArgs(valid_args=valid_args).arg_dict
        assert len(arg_dict) > 0, usage

        point_list = None

        # load the master processor
        masterProc = master.LstarMaster(comm=comm)

        if arg_dict['filename'] is not None:
            # load the filename to the master process
            # TODO: Need to adjust this to work on a list of files
            masterProc.add_file_time_pair(tuple((arg_dict['filename'], 0, arg_dict['vector'])))
        else:
            # We need a file
            print("Cannot continue with a filename")
            assert False

        if arg_dict['radius'] is not None:
            # get list of points within the RE distance requested
            # TODO: Need to work on list of file names
            data = gpd.VtkData(filename=arg_dict['filename'], vector=arg_dict['vector'])
            points = np.array(data.points)
            bRE = float(arg_dict['radius'])
            pRE = algc.mag(points)
            pREu = pRE < bRE
            pREl = pRE > data.get_calc_boundary()
            pREsel = np.where(pREu == pREl)
            point_list = points[pREsel]

            master_list = []
            for p in point_list:
                pdata = np.concatenate([p, [0]])
                master_list.append(pdata)

        elif arg_dict['pointsFile'] is not None:
            # get list of points/times from file
            print ("Getting Points list from file...")
            assert False, "File Reading not yet implemented."
        else:
            print ("Must have a list of points to process, or a radius to process")
            assert False

        # load points to the master process
        for pair in master_list:
            masterProc.add_point_time_pair(pair)

        print ("Number of points to process: {}".format(masterProc.get_number_of_points()))

        # start master process
        masterProc()

    else:
        ####################
        #   WORKER RANKS   #
        ####################
        # Worker dependent modules
        import worker
        workerProc = worker.LstarWorker(comm=comm)
        workerProc()


if __name__ == '__main__':
    sys.exit(main())