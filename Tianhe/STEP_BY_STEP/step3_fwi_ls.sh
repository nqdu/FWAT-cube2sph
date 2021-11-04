#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:59:59
#SBATCH --job-name LINESEARCH
#SBATCH --output=LINESEARCH_%j.txt
#SBATCH --partition=debug
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
#SBATCH --mail-type=ALL

set -e 

#####################################
###### Parameters ##################

# event sets
setb=12
sete=22

# simu_type
simu_type=tele 

# start model
startmod=M00

# iterations 
iters=1

# HPC facilities
nnodes=4
ntasks_per_node=40

# FKSEM path
fksem='/home/l/liuqy/nqdu/specfem3d/'

###### Parameters END ################
####################################

# LOAD all module
source activate base
module load NiaEnv/2019b
module load intel openmpi 

echo "========================="
echo "Line Search begin `date`"
echo "========================="

# machinefile 
NPROC=`grep "^NPROC" DATA/Par_file | awk '{print $3}'`
srun hostname | sort -n | uniq | awk -F. '{print $1}' > slurm.host 
python utils/generate_hostfile.py $nnodes $ntasks_per_node $NPROC
nodes=`ls slurm.host.? | wc -l`

# go with tele data 
startidx=$(printf %d `echo $startmod | cut -d'M' -f2`)
for((ii=0;ii<$iters;ii++));
do  
    # change LOCAL_PATH to the model database
    jj=$(printf %02d $(echo "$ii+$startidx"|bc))

    # linesearch
    echo "Line Searching ..."
    if [[ $simu_type == "tele" ]];then 
        cp src_rec/sources_ls.dat.tele src_rec/sources_ls -r 
    else
        cp src_rec/sources_ls.dat.noise src_rec/sources_ls -r 
    fi 
    python src/run_linesearch.py M$jj $simu_type $NPROC $nodes

done 

echo "========================="
echo "Line Search finished `date`"
echo "========================="
