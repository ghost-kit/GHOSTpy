from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

data = (rank+1)**2
data = comm.gather(data, root=0)
if rank == 0:
    for i in range(size):
        assert data[i] == (i+1)**2
        print (i,  ": True")

        if i != 0:
            comm.send("Test Data {}".format(i), dest=i)

else:
    assert data is None
    test = comm.recv(source=0)
    print ("Test Received: {}".format(test))