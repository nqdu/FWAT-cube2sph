import os 
import sys 

def AllocateTasks(ntasks_per_node,nnodes,nproc_per_evt,nprocs,lines):
    ntasks_remain = [ntasks_per_node for i in range(nnodes)]
    hostname = [[] for i in range(nprocs)]

    for i in range(nprocs):
        n = nproc_per_evt
        for j in range(nnodes):
            if ntasks_remain[j] ==0:
                continue 
            if n > ntasks_remain[j]:
                n -= ntasks_remain[j]
                hostname[i].append(f"{lines[j].split()[0]} slots={ntasks_remain[j]}\n")
                ntasks_remain[j] = 0
            elif n == ntasks_per_node:
                ntasks_remain[j] = 0
                hostname[i].append(f"{lines[j].split()[0]} slots={n}\n")
                n = 0
            else:
                ntasks_remain[j] -= n 
                hostname[i].append(f"{lines[j].split()[0]} slots={n}\n")
                n = 0
            if n==0:
                break 
    
    return hostname

def AllocateTasksBalance(ntasks_per_node,nnodes,nproc_per_evt,nprocs,lines):
    ntasks_remain = [ntasks_per_node for i in range(nnodes)]
    hostname = [[] for i in range(nprocs)]

    # compute minimum 


def main():
    if len(sys.argv) !=4:
        print("usage: python thisfile nnodes ntasks_per_node,nprocs_per_evt")
        exit(1)
    
    # get input args
    nnodes = int(sys.argv[1])
    ntasks_per_node = int(sys.argv[2])
    nproc_per_evt = int(sys.argv[3])
    
    # read hostfile
    f = open("slurm.host","r")
    lines = f.readlines()
    f.close()
    assert(len(lines) == nnodes)

    # compute # of simultaneous tasks
    nprocs = (nnodes * ntasks_per_node) // nproc_per_evt
    if nprocs ==0 :
        print("please increase your nodes!")
        exit(1)
    else:
        print("# of simultaneous tasks = %d" %nprocs)

    # allocate tasks
    hostname = AllocateTasks(ntasks_per_node,nnodes,nproc_per_evt,nprocs,lines)

    # write to slurm.new
    for i in range(nprocs):
        f = open(f"slurm.host.{i}","w")
        f.writelines(hostname[i])
        f.close()
    print(hostname)

main()
