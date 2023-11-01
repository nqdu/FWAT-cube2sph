#!/bin/bash
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=24
#SBATCH --ntasks-per-node=40
#SBATCH --time=02:02:59
#SBATCH --job-name LS
#SBATCH --output=LS%j_%a.txt
#SBATCH --mem=0

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
#module load NiaEnv/2019b
module load intel openmpi

. utils/parameters.sh 

NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`
sizes=(`grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`)
model=M20
MODEL=${model}_step${sizes[$SLURM_ARRAY_TASK_ID]}
SET=ls
SIMU_TYPE=rf


\cp DATA/Par_file.$SIMU_TYPE DATA/Par_file
\cp fwat_params/FWAT.PAR.$SIMU_TYPE fwat_params/FWAT.PAR
\cp fwat_params/MEASUREMENT.PAR.$SIMU_TYPE fwat_params/MEASUREMENT.PAR
\cp src_rec/sources_ls.dat.$SIMU_TYPE src_rec/sources_ls.dat -r 

# substitute
LOCAL_PATH="./OUTPUT_FILES/DATABASES_MPI"
./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH DATA/Par_file
./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH  DATA/meshfem3D_files/Mesh_Par_file
mpirun -np $NPROC $fksem/bin/xfwat3_linesearch $MODEL $SET $SIMU_TYPE

mv output_fwat3_log_${model}_step${sizes[$SLURM_ARRAY_TASK_ID]}.ls.txt output_fwat3_log_${model}_${SIMU_TYPE}_step${sizes[$SLURM_ARRAY_TASK_ID]}.ls.txt
