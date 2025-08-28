import yaml 
from .system.system import sys_remove
from .system.specfem import get_param

def run(args:list[str]):
    """
    clean the simulation directory for a given model and event id

    Parameters
    ----------
    args : list[str]
        list of arguments, should contain two or three elements:
        1. model name (e.g., M00)
        2. event id (e.g., P0001)
        3. (optional) deepclean (True/False), if True, remove the entire model directory 
    """
    if len(args) !=3 and len(args) !=2 :
        print("Usage: fwat clean MODEL(M00) evtid (deepclean)")
        exit(1)

    # get args
    mod = args[0]
    evtid = args[1]
    deepclean = False 
    if len(args) == 3:
        deepclean = bool(args[2])
    
    # get solver location
    from fwat.const import SOLVER
    rundir = f"{SOLVER}/{mod}/{evtid}"

    # get local path
    LOCAL_PATH = get_param(f"{rundir}/DATA/Par_file","LOCAL_PATH")
    
    # remove
    sys_remove(f"{rundir}/{LOCAL_PATH}")
    sys_remove(f"{rundir}/OUTPUT_FILES/*.sac")
    sys_remove(f"{rundir}/DATA/meshfem*")
    sys_remove(f"{rundir}/OUTPUT_FILES/*forward_wavefield.bin*")

    # remove period band
    sys_remove(f"{rundir}/OUTPUT_FILES/T[0-9]*_T[0-9]*")

    # remove timestamp
    sys_remove(f"{rundir}/OUTPUT_FILES/timestamp*")

    # deep clean
    if deepclean:
        sys_remove(f"{rundir}")