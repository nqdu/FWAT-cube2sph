from multiprocessing import Pool
import sys 
import os 
import numpy as np
from parameters import *


def run_fwdadj(arg):
    rank,setnum,nproc_per_src,mod,simu_type= arg 

    cmd=f"mpirun -machinefile slurm.host.{rank} --oversubscribe -np {nproc_per_src} "\
          f"{fksem}/bin/xfwat0_forward_data {mod} set{setnum} {simu_type} > forward.set{setnum}.txt"
    print(cmd,flush=True)
    os.system(cmd)

def main():
    if len(sys.argv) != 7:
        print("usage python forward.py mod setb sete simu_type nproc nodes")
        print("example:")
        print("python forward.py M00 2 20 tele 80 4")
        exit(1)
    else:
        mod= sys.argv[1]
        setb = int(sys.argv[2])
        sete = int(sys.argv[3])
        simu_type = sys.argv[4]
        nproc_per_src = int(sys.argv[5])
        nodes = int(sys.argv[6])
    
    # mkdir 
    os.system("mkdir -p data solver")

    # arglist
    arglist = []
    rank = 0
    for i in range(setb,sete+1):
        arglist.append([rank,i,nproc_per_src,mod,simu_type])
        rank += 1
        if rank == nodes or i==sete: # all tasks are in queue
            # run 
            pool = Pool(len(arglist))
            pool.map(run_fwdadj,arglist)
            pool.close()
            pool.join()
            
            # renew arglist and rank
            rank = 0
            arglist = []


main()