#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:59:59
#SBATCH --job-name POST
#SBATCH --output=POST_%j.txt
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
echo "POST-PROCESSING begin `date`"
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
    jj=$(printf %02d ` echo "$ii + $startidx" | bc`)
    mkdir -p optimize
    if [ $jj == "00" ]; then 
        rm -rf optimize/MODEL_M00
        cp -r initial_model optimize/MODEL_M00
    fi 
    echo " mpirun -oversubscribe -np $NPROC $fksem/bin/xfwat2_postproc_opt M$jj set$setb set$sete true "
    mpirun -oversubscribe -np $NPROC $fksem/bin/xfwat2_postproc_opt M$jj set$setb set$sete true > POST 
    for step in `grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`;do
        echo "======================================="
        echo  Meshing for model $MODEL step:$step
        echo "========================================"
        #=====
        ./utils/change_par_file.sh LOCAL_PATH ./optimize/MODEL_M${jj}_step${step} DATA/Par_file
        ./utils/change_par_file.sh LOCAL_PATH ./optimize/MODEL_M${jj}_step${step} DATA/meshfem3D_files/Mesh_Par_file
        ./utils/change_par_file.sh SAVE_MESH_FILES .false. DATA/Par_file       
        mpirun -machinefile slurm.host -oversubscribe -np $NPROC $fksem/bin/xmeshfem3D
        mpirun -machinefile slurm.host -oversubscribe -np $NPROC $fksem/bin/xgenerate_databases
    done

done 

echo "========================="
echo "POST-PROCESSING finished `date`"
echo "========================="
