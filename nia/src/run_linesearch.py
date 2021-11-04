from multiprocessing import Pool
import sys 
import os 
from parameters import *

def run_ls(arg):
    rank,step,nproc_per_src,mod,simu_type = arg 

    cmd=f"mpirun -machinefile slurm.host.{rank} --oversubscribe -np {nproc_per_src}" \
        f" {fksem}/bin/xfwat3_linesearch {mod}_step{step} ls {simu_type} > LS.{step}.txt "
    print(cmd,flush=True)
    os.system(cmd)

def main():
    if len(sys.argv) != 5:
        print("usage python run_ls.py mod simu_type nproc nodes")
        exit(1)
    else:
        mod = sys.argv[1]
        simu_type = sys.argv[2]
        nproc_per_src = int(sys.argv[3])
        nodes = int(sys.argv[4])
    
    # mkdir 
    os.system("mkdir -p solver misfits")

    # get line search parameters
    f = os.popen("grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'")
    line = f.readline()
    f.close()
    steps = line.split()

    # arglist
    arglist = []
    rank = 0
    for i in range(len(steps)):
        arglist.append([rank,steps[i],nproc_per_src,mod,simu_type])
        rank += 1
        if rank == nodes or i== len(steps)-1:
            # run
            pool = Pool(len(arglist))
            pool.map(run_ls,arglist)
            pool.close()
            pool.join()
            
            # renew arglist and rank
            rank = 0
            arglist = []
main()