#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:59:59
#SBATCH --job-name FWDADJ
#SBATCH --output=FWDADJ_%j.txt
#SBATCH --partition=TH_HPC3

set -e 

#####################################
###### Parameters ##################

# event sets
setb=12
sete=22

# simu_type
simu_type=tele 

# start model
startmod=M01

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
echo "Forward and Adjoint simulation begin `date`"
echo "========================="

# copy parameters if required
if [ $simu_type == 'noise' ];then
    echo "copying params for noise"
    cp DATA/Par_file.noise DATA/Par_file
    cp fwat_params/FWAT.PAR.noise fwat_params/FWAT.PAR
    cp fwat_params/MEASUREMENT.PAR.noise fwat_params/MEASUREMENT.PAR
    cp src_rec/sources_ls.dat.noise src_rec/sources_ls.dat
elif [ $simu_type == 'tele' ];then
    echo "copying params for tele"
    cp DATA/Par_file.tele DATA/Par_file
    cp fwat_params/FWAT.PAR.tele fwat_params/FWAT.PAR
    cp fwat_params/MEASUREMENT.PAR.tele fwat_params/MEASUREMENT.PAR
    cp src_rec/sources_ls.dat.tele src_rec/sources_ls.dat
fi 

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
    jj=$(printf %02d ` echo "$ii + $startidx" | bc`)
    if [ $jj == "00" ];then
        ./utils/change_par_file.sh LOCAL_PATH ./OUTPUT_FILES/DATABASES_MPI DATA/Par_file
        ./utils/change_par_file.sh LOCAL_PATH ./OUTPUT_FILES/DATABASES_MPI DATA/meshfem3D_files/Mesh_Par_file
    else
        ./utils/change_par_file.sh LOCAL_PATH ./optimize/MODEL_M${jj} DATA/Par_file
        ./utils/change_par_file.sh LOCAL_PATH ./optimize/MODEL_M${jj} DATA/meshfem3D_files/Mesh_Par_file
    fi 
    
    # forward and adjoint simulation
    echo "Forward and Adjoint simulations ..."
    python src/run_fwdadj.py M$jj $setb $sete $simu_type $NPROC $nodes 

done 

echo "========================="
echo "Forward and Adjoint simulation finished `date`"
echo "========================="
