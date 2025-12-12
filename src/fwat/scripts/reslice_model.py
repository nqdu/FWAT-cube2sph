from mpi4py import MPI
import numpy as np
from numba import jit  
import glob 

from fwat.FortranIO import FortranIO

@jit(nopython=True)
def coords2discon(xstore,ystore,zstore,ibool):
    npts = len(ibool)
    x = np.zeros(npts,dtype=np.int32)
    y = np.zeros(npts,dtype=np.int32)
    z = np.zeros(npts,dtype=np.int32)

    for i in range(npts):
        idx = ibool[i]
        x[i] = xstore[idx]
        y[i] = ystore[idx]
        z[i] = zstore[idx]
    return x,y,z

def read_coordinates(filename):
    f = FortranIO(filename)
    _ = f.read_record('i4')[0] # nspec
    _ = f.read_record('i4')[0] #gnlob
    _ = f.read_record('i4')
    ibool = f.read_record('i4') - 1
    xstore = f.read_record('f4')
    ystore = f.read_record('f4')
    zstore = f.read_record('f4')

    x,y,z = coords2discon(xstore,ystore,zstore,ibool)

    return x,y,z 

def run(argv):
    if len(argv) != 3:
        print("Usage: fwat-utils reslice input_dir output_dir param")
        print("\texample: mpirun -np 100 fwat-utils reslice DATABASES_MPI_100cores DATABASES_MPI_4cores vp")
        exit(1) 
    
    # init MPI 
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # get args 
    input_dir = argv[0]
    output_dir = argv[1]
    param = argv[2]
    if myrank ==0:
        print(f"input_dir = {input_dir}, output_dir = {output_dir}, param = {param}")

    # get info from input file 
    filename = input_dir + '/proc%06d'%(myrank) + '_external_mesh.bin'
    xin,yin,zin = read_coordinates(filename)

    # read model file 
    filename = input_dir + '/proc%06d'%(myrank) + f'_{param}.bin'
    f = FortranIO(filename)
    model_in = f.read_record('f4')
    f.close()

    # create dictionary for model_in 
    model_dict = {}
    for i in range(len(xin)):
        model_dict[(xin[i],yin[i],zin[i])] = model_in[i]
    
    # get how many files in output_dir 
    nprocs_out = len(glob.glob(f"{output_dir}/proc*_external_mesh.bin"))
    if myrank == 0:
        print(f"number of files to interp = {nprocs_out}")
    
    # loop each file to interp
    for iproc_out in range(nprocs_out):
        if myrank ==0:
            print(f"interp for file {iproc_out} ...")
        
        # get info 
        filename = output_dir + '/proc%06d'%(iproc_out) + '_external_mesh.bin'
        xout,yout,zout = read_coordinates(filename)
        if myrank ==0:
            print(f"finish read coordinates ...")

        # loop each point in output to find nearest in input
        dist_idx = np.zeros((len(xout)),dtype=np.int32)
        for i in range(len(xout)):
            key = (xout[i],yout[i],zout[i])
            if key in model_dict:
                dist_idx[i] = i 
            else:
                dist_idx[i] = -1
        comm.barrier()


        if myrank ==0:
            print(f"trying to find closest point ...")

        model_out = np.zeros((len(xout)),dtype=np.float64)
        count = np.zeros((len(xout)),dtype=np.int32)
        for i in range(len(xout)):
            if dist_idx[i] != -1:
                model_out[i] = model_dict[(xout[i],yout[i],zout[i])]
                count[i] += 1
            else:
                model_out[i] = 0.0
        model_out = comm.allreduce(model_out,op=MPI.SUM)
        count = comm.allreduce(count,op=MPI.SUM) 

        idx = np.where(count > 0) 
        model_out[idx] = model_out[idx] / count[idx]

        # save output 
        filename = output_dir + '/proc%06d'%(iproc_out) + f'_{param}.bin'
        if myrank == 0: 
            print(f"writing to file {filename} ")
            print("=======================\n")

        model_out = np.asarray(model_out,dtype='f4')
        f = FortranIO(filename,'w') 
        f.write_record(model_out)
        f.close()


def main():
    import sys

    run(sys.argv[1:])


if __name__ == "__main__":
    main()
