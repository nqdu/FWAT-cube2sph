import numpy as np
import os
from string import Template
from mpi4py import MPI
import h5py

def _get_first_dim_npy(filename):
    """Read the first dimension of a numpy binary file (.npy).
    This function reads the header of a numpy binary file to determine the size of the first dimension.
    It is useful for files that are not too large to fit into memory.

    Parameters
    ----------
    filename : str
        The path to the numpy binary file.

    Returns
    -------
    int
        The size of the first dimension of the array stored in the file.
"""
    import numpy.lib.format

    with open(filename,"rb") as fobj:
        version = numpy.lib.format.read_magic(fobj)
        func_name = 'read_array_header_' + '_'.join(str(v) for v in version)
        func = getattr(numpy.lib.format, func_name)
        
        # get header
        out = func(fobj)

    return out[0][0]


def _count_lines(filename):
    """
    Count the number of lines in a file, handling the case where the last line may not end with a newline character.
    
    Parameters
    ----------
    filename : str
        The path to the file.

    Returns
    -------
    int
        The number of lines in the file.
    
    """
    with open(filename, 'rb') as f:
        count = 0
        while True:
            buf = f.read(1024 * 1024)
            if not buf:
                break
            count += buf.count(b'\n')
            last_buf = buf
        # Check if the last line is missing a newline
        if last_buf and not last_buf.endswith(b'\n'):
            count += 1
        return count
    
def _read_rotation_file(filename):
    """
    Read a rotation matrix file and return a dictionary with keys as station names and values as rotation matrices.
    The file format is expected to have the first line with network and station name (e.g., TA A01) followed by three lines representing a 3x3 matrix.
    
    Parameters
    ----------
    filename : str
        The path to the rotation matrix file.

    Returns
    -------
    dict
        A dictionary where keys are station names (e.g., "TA.A01") and values are 3x3 numpy arrays representing the rotation matrices.
    """
    data = {}
    with open(filename) as f:
        lines = [line.strip() for line in f if line.strip()]  # remove blank lines

    i = 0
    while i < len(lines):
        # First line contains two strings, e.g., TA     A01
        parts = lines[i].split()
        key = f"{parts[0]}.{parts[1]}"  # e.g., TA.A01
        i += 1

        # Next 3 lines are the 3x3 matrix
        matrix = []
        for _ in range(3):
            row = list(map(float, lines[i].split()))
            matrix.append(row)
            i += 1

        data[key] = np.array(matrix, dtype=float)

    return data


def rotate_seismo_fwd(
        rot_list:list,
        fn_matrix:str,
        from_dir:str,to_dir:str,
        from_template_str:str,
        to_template_str:str):
    """
    Function to rotate the forward seismograms from the Cube2sph system to the local Geographic system, MPI parallelized.

    Parameters
    ----------
    fn_matrix : str
        Path to the rotation matrix file produced by the write_station_file program. Rows: (N, E, Z) in geographic system; Columns: (x, y, z) in Cartesian system.
    from_dir : str
        Input directory containing the seismograms.h5 in the Cube2sph system.
    to_dir : str
        Output directory where the rotated seismograms will be saved in the local Geographic system, saved as .npy files.
    from_template_str : str
        File name template for input seismograms. Should include placeholders for network, station, and component, e.g., '${nt}.${sta}.BX${comp}.semd'.
    to_template_str : str
        File name template for output seismograms. Should include placeholders for network, station, and component, e.g., '${nt}.${sta}.BX${comp}.sem.ascii'.

    Returns
    -------
    None
    """
    # input is hdf5
    #--MPI-
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nproc = comm.Get_size()

    from_template = Template(from_template_str)
    to_template = Template(to_template_str)

    rotate = "XYZ->NEZ"
    comp_left = rotate[0:3]
    comp_right = rotate[5:8]
    from_comp, to_comp = comp_left, comp_right

    # read fn matrix
    rot_dict = _read_rotation_file(fn_matrix)

    # open hdf5
    fio = h5py.File(from_dir + "/seismograms.h5","r")

    # get station names
    keys = rot_list
    nrec_loc = len(keys)

    # allocate jobs
    seismo_dict = {}

    for i_sta in range(nrec_loc):
        nt,sta = keys[i_sta].split('.')
        nu = rot_dict[keys[i_sta]]
        missing_file = False
        nstep = -1
        for i_comp in range(0, 3):
            if (from_comp[i_comp] == '0'):
                continue 
            dname = from_template.substitute(nt=nt, sta=sta, comp=from_comp[i_comp])
            if dname not in fio.keys():
                #print(f"{dname} does not exist but required by rotation, skipping this station")
                missing_file = True
                break

            if nstep < 0:
                nstep = fio[dname].shape[0]
        
        if (missing_file):
            continue
        seis = np.zeros(shape=(nstep, 3), dtype=float)
        arr = np.zeros(shape=(nstep, 2), dtype=float)
        for i_comp in range(0, 3):
            if (from_comp[i_comp] == '0'):
                continue

            dname = from_template.substitute(nt=nt, sta=sta, comp=from_comp[i_comp])
            arr[:] =  np.asarray(fio[dname][:])
            seis[:,i_comp] = arr[:,1]


        # apply rotation matrix
        seis = np.matmul(seis, np.transpose(nu))
        for i_comp in range(0, 3):
            if (to_comp[i_comp] == '0'):
                continue
            fn = os.path.join(to_dir, 
                    to_template.substitute(nt=nt, sta=sta, comp=to_comp[i_comp]))
            #print(fn)
            arr[:,1] = seis[:,i_comp]
            #print(f"writing to ${fn}")
            seismo_dict[fn] = arr * 1.

    comm.barrier()
    fio.close()

    return seismo_dict

def rotate_seismo_adj(
        rot_list:list,
        fn_matrix:str,from_dir:str,to_dir:str,
        from_template_str:str,to_template_str:str):
    """
    Function to rotate the adjoint sources from the local Geographic system back to the Cube2sph system, MPI parallelized.
    
    Parameters
    ----------
    fn_matrix : str
        Path to the rotation matrix file produced by the write_station_file program.
    from_dir : str
        Input directory containing the *.adj.npy in the local Geographic system.
    to_dir : str
        Output directory where the rotated seismograms will be saved in the Cube2sph system, saved as .txt file.
    from_template_str : str
        File name template for input seismograms. Should include placeholders for network, station, and component, e.g., '${nt}.${sta}.BX${comp}.sem.ascii'.
    to_template_str : str
        File name template for output seismograms. Should include placeholders for network, station, and component, e.g., '${nt}.${sta}.BX${comp}.semd'.
    
    """

    #--MPI-
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nproc = comm.Get_size()

    from_template = Template(from_template_str)
    to_template = Template(to_template_str)

    rotate = "XYZ<-NEZ"
    comp_left = rotate[0:3] 
    comp_right = rotate[5:8]
    from_comp, to_comp = comp_right, comp_left

    #arr = np.zeros(shape(nsteps, 2), dtype=float)
    # read fn matrix
    rot_dict = _read_rotation_file(fn_matrix)
    keys = rot_list
    nrec_local = len(keys)

    for i_sta in range(0, nrec_local):
        nt,sta = keys[i_sta].split('.')
        nu = rot_dict[keys[i_sta]]
        missing_file = False
        nstep = -1
        for i_comp in range(0, 3):
            if (from_comp[i_comp] == '0'):
                continue 
            dname = from_template.substitute(nt=nt, sta=sta, comp=from_comp[i_comp])
            fn = os.path.join(from_dir, 
                    from_template.substitute(nt=nt, sta=sta, comp=from_comp[i_comp]))
            if (not os.path.isfile(fn)):
                #print(fn)
                #print(f"{dname} does not exist but required by rotation, skipping this station")
                #os.system(f"wc -l {fn}")
                missing_file = True
                break
            if nstep < 0:
                nstep = _get_first_dim_npy(fn)
        
        if (missing_file):
            continue
        seis = np.zeros(shape=(nstep, 3), dtype=float)
        arr = np.zeros(shape=(nstep, 2), dtype=float)
        for i_comp in range(0, 3):
            if (from_comp[i_comp] == '0'):
                continue
            fn = os.path.join(from_dir, 
                from_template.substitute(nt=nt, sta=sta, comp=from_comp[i_comp]))
            arr = np.load(fn)
            seis[:,i_comp] = arr[:,1]
        
        seis = np.matmul(seis, nu)

        for i_comp in range(0, 3):
            if (to_comp[i_comp] == '0'):
                continue
            fn = os.path.join(to_dir, 
                    to_template.substitute(nt=nt, sta=sta, comp=to_comp[i_comp]))
            #print(fn)
            arr[:,1] = seis[:,i_comp]
            #print(f"writing to ${fn}")
            #arr.tofile(fn)
            np.savetxt(fname=fn, X=arr, fmt='%11.6f%19.7E')

    comm.barrier()