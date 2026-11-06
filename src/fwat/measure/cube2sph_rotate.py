import numpy as np
import os
from string import Template
from mpi4py import MPI
import h5py

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
        Namespace prefix for keys of the rotated traces returned in memory.
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
        fn_matrix:str,seismo_adj:dict[str,np.ndarray],
        from_dir:str,to_dir:str,
        from_template_str:str,to_template_str:str):
    """
    Function to rotate the adjoint sources from the local Geographic system back to the Cube2sph system, MPI parallelized.
    
    Parameters
    ----------
    fn_matrix : str
        Path to the rotation matrix file produced by the write_station_file program.
    seismo_adj : dict[str, np.ndarray]
        In-memory local geographic adjoint traces, keyed by their source paths.
        Each trace has time and amplitude columns.
    from_dir : str
        Directory used to form keys in seismo_adj.
    to_dir : str
        Output directory where the rotated seismograms will be saved in the Cube2sph system, saved as .txt file.
    from_template_str : str
        Dictionary key template for input traces, e.g., '${nt}.${sta}.BX${comp}.adj.sem.npy'.
    to_template_str : str
        File name template for solver output, e.g., '${nt}.${sta}.BX${comp}.adj'.
    
    """

    #--MPI-
    comm = MPI.COMM_WORLD

    from_template = Template(from_template_str)
    to_template = Template(to_template_str)

    rotate = "XYZ<-NEZ"
    comp_left = rotate[0:3] 
    comp_right = rotate[5:8]
    from_comp, to_comp = comp_right, comp_left

    # read fn matrix
    rot_dict = _read_rotation_file(fn_matrix)
    keys = rot_list
    nrec_local = len(keys)

    for i_sta in range(0, nrec_local):
        nt,sta = keys[i_sta].split('.')
        nu = rot_dict[keys[i_sta]]
        input_names = []
        for i_comp in range(0, 3):
            if (from_comp[i_comp] == '0'):
                continue 
            input_names.append(os.path.join(
                from_dir,
                from_template.substitute(nt=nt, sta=sta, comp=from_comp[i_comp])))

        if any(fn not in seismo_adj for fn in input_names):
            continue
        nstep = seismo_adj[input_names[0]].shape[0]
        seis = np.zeros(shape=(nstep, 3), dtype=float)
        for i_comp, fn in enumerate(input_names):
            arr = seismo_adj[fn]
            if arr.shape != (nstep, 2):
                raise ValueError(f"Invalid adjoint trace shape for {fn}: {arr.shape}")
            seis[:,i_comp] = arr[:,1]
        
        seis = np.matmul(seis, nu)

        # Use the time column from the last input, as in the file-based path.
        arr = seismo_adj[input_names[-1]].copy()

        for i_comp in range(0, 3):
            if (to_comp[i_comp] == '0'):
                continue
            fn = os.path.join(to_dir, 
                    to_template.substitute(nt=nt, sta=sta, comp=to_comp[i_comp]))
            arr[:,1] = seis[:,i_comp]
            with open(fn, 'w') as fio:
                fio.write('\n'.join([f"{row[0]:11.6f} {row[1]:19.7E}" for row in arr]))


    comm.barrier()
