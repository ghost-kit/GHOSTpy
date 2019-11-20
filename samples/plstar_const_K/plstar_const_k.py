from mpi4py import MPI
import sys
import lsArgs
from ghostpy.data import VtkData as vtkd
from ghostpy.data import GpData as gpd

import worker
from ghostpy.algorithms import common as algc
import numpy as np

# setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def main():
    if rank == 0:
        print ("Output is being re-directed to rank-dependent files.  Please see these files for output.")

    sys.stdout = open("stdout_rank{}.txt".format(rank), "w+", buffering=1)
    sys.stderr = open("stderr_rank{}.txt".format(rank), "w+", buffering=1)

    # get the arguments from the command line
    usage = "System Usage: {} -f <filename> -v <vector> -p <pointsFile> -r <pointsRadius> -t <convergence tolerance>".format(sys.argv[0])
    valid_args = {'f': "filename", 'v': "vector", 'r': "radius", 'p': "pointsFile", 't': "tol"}
    arg_dict = lsArgs.lsArgs(valid_args=valid_args).arg_dict
    assert len(arg_dict) > 0, usage

    tol = float(arg_dict['tol'])

    worker_responsibility = {}
    data = vtkd.VtkData(filename=arg_dict['filename'], vector=arg_dict['vector'])
    points = np.array(data.points)
    bRE = float(arg_dict['radius'])
    pRE = algc.mag(points)
    pREu = pRE < bRE
    pREl = pRE > data.get_calc_boundary()
    pREsel = np.where(pREu == pREl)
    point_list = points[pREsel]

    parts = len(point_list)

    num_workers = size
    points_per_proc = float(parts) / num_workers

    if rank == 0:
        print ("Total Number of Points: {}".format(parts))
        print("Maximum number of points per processor: {}".format(np.ceil(points_per_proc)))
        print("Convergence Tolerance: {}".format(tol))

    for proc in range(num_workers):
        worker_responsibility[proc] = []

    # parting out points

    for idx in range(parts):
        procNum = idx % num_workers
        worker_responsibility[procNum].append(point_list[idx])

    workerbot = worker.LStarConvergeWorker(comm=comm)

    workerbot(conv_tol=tol, data=data, vector=arg_dict['vector'], pointlist=worker_responsibility[rank], k=0)

    data = comm.gather(workerbot.get_data(), root=0)

    if rank == 0:
        grid_dat = vtkd.VtkData(filename=arg_dict['filename'], vector='B')
        lu = ReverseLookup(data=grid_dat)

        locx = np.array([])
        locy = np.array([])
        locz = np.array([])
        lstar = np.array([])
        clstar = np.array([])
        clines = np.array([])

        for datum in data:
            locx = np.concatenate([locx, (datum[0])[:, 0]])
            locy = np.concatenate([locy, (datum[0])[:, 1]])
            locz = np.concatenate([locz, (datum[0])[:, 2]])
            lstar = np.concatenate([lstar, datum[1]])
            clstar = np.concatenate([clstar, datum[2]])
            clines = np.concatenate([clines, datum[3]])

        lines = __format_data_for_output__(filename=arg_dict['filename'],
                                           vector=arg_dict['vector'],
                                           locx=locx,
                                           locy=locy,
                                           locz=locz,
                                           lstar=lstar,
                                           clstar=clstar,
                                           clines=clines,
                                           tol=tol)

        lstar_aligned = np.zeros(shape=[len(grid_dat.points)])
        lstar_aligned[:] = np.NaN

        lstarc_aligned = np.zeros(shape=[len(grid_dat.points)])
        lstarc_aligned[:] = np.NaN

        clines_aligned = np.zeros(shape=[len(grid_dat.points)])
        clines_aligned[:] = np.NaN

        for idx in range(len(locx)):
            pt = tuple([locx[idx], locy[idx], locz[idx]])
            lstar_aligned[lu(pt)] = lstar[idx]
            lstarc_aligned[lu(pt)] = clstar[idx]
            clines_aligned[lu(pt)] = clines[idx]

        # print (lstarc_aligned)
        # for val in lstarc_aligned:
        #     print (val)

        print (len(np.argwhere(np.isfinite(lstarc_aligned))))

        mdims = grid_dat.dims[:-1].copy()
        mdims[0] = grid_dat.dims[2]
        mdims[2] = grid_dat.dims[0]

        gx,gy,gz = build_grid(grid_dat.points, mdims)

        ldat = build_dat(lstar_aligned, mdims)
        lcdat = build_dat(lstarc_aligned, mdims)
        llines = build_dat(clines_aligned, mdims)

        from pyevtk.hl import gridToVTK
        gridToVTK("lfm_lstar_test", gx, gy, gz, pointData={"L": ldat, "LC": lcdat, "Lines": llines})

        ofile = open("outdata.csv", "w+")
        ofile.writelines(lines)


def __format_data_for_output__(filename, vector, locx, locy, locz, lstar, clstar, clines, tol):
    lines = ["File Name,\t vector,\t SM Loc (X), \t SM Loc (Y), \t SM Loc (Z), "
             "\t L*,\t converged L*,\t Convergence Lines \t, tol\n"]

    for idx in range(len(locx)):
        dLine = str(filename) + ",\t" + str(vector) + ",\t" + str(locx[idx]) + ",\t" + str(locy[idx]) + ",\t" + \
                str(locz[idx]) + ",\t" + str(lstar[idx]) + ",\t" + str(clstar[idx]) + ",\t" + \
                str(clines[idx]) + ",\t" + str(tol) + "\n"

        lines.append(dLine)

    # for dat in input_dict:
    #     filename = self.fileTimeList[np.array(dat)[-1]][0]
    #     vector = self.fileTimeList[np.array(dat)[-1]][1]
    #     ddat = np.array(dat)[:-1]
    #     ddat = str(ddat)
    #     ddat = ddat.replace("[", "")
    #     ddat = ddat.replace("]", "")
    #     ddat = re.sub("\s+", ",\t ", ddat.strip())
    #     d = input_dict[dat]
    #     dline = str(filename) + ",\t" + str(vector) + ",\t" + str(ddat) + ",\t" + str(d) + '\n'
    #     lines.append(dline)

    return lines


def build_grid(data, dims):
    xyz = np.array(data)
    x,y,z = np.hsplit(xyz,3)
    x = x.flatten()
    y = y.flatten()
    z = z.flatten()

    x = x.reshape(dims)
    y = y.reshape(dims)
    z = z.reshape(dims)

    return x, y, z


def build_dat(dat, dims):
    return dat.reshape(dims)


class ReverseLookup(object):
    def __init__(self, data):
        assert isinstance(data, vtkd.VtkData)
        points = data.points
        self.lu = {}

        for idx in range(len(points)):
            self.lu[tuple(points[idx])] = idx

    def __call__(self, point):
        return self.lu[tuple(point)]






if __name__ == '__main__':
    sys.exit(main())


