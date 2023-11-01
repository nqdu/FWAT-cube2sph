#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --array=21-26
#SBATCH --time=01:48:59
#SBATCH --job-name FWD_ADJ
#SBATCH --output=FWD_ADJ-%j_set%a.txt
#SBATCH --mem=0

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
#module load NiaEnv/2019b
module load intel openmpi

# include file
. utils/parameters.sh

#==== Comment out the following if running SEM mesh with new models====#
MODEL=M20
simu_type=rf

# copy 
echo "copying params for $simu_type"
\cp DATA/Par_file.$simu_type DATA/Par_file
\cp fwat_params/FWAT.PAR.$simu_type fwat_params/FWAT.PAR
\cp fwat_params/MEASUREMENT.PAR.$simu_type fwat_params/MEASUREMENT.PAR

mod=$MODEL
if [ ${mod} == "M00"  ]; then
  LOCAL_PATH="./OUTPUT_FILES/DATABASES_MPI"
  ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH DATA/Par_file
  ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH  DATA/meshfem3D_files/Mesh_Par_file
else
  LOCAL_PATH="./optimize/MODEL_${mod}"
  ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH DATA/Par_file
  ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH  DATA/meshfem3D_files/Mesh_Par_file
fi

# SET
SET=set$SLURM_ARRAY_TASK_ID
set_num=`echo $SET |awk -Fset '{printf"%d\n",$2}'`
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`

# set_num > 10 for tele
mpirun -np $NPROC $fksem/bin/xfwat1_fwd_measure_adj $MODEL $SET $simu_type 3
