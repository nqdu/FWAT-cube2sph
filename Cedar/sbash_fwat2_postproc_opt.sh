#!/bin/bash
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:55:59
#SBATCH --job-name=POST
#SBATCH --output=POST_%j.txt
#SBATCH --mem=0

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
#module load NiaEnv/2019b
module load intel openmpi

# include file
. utils/parameters.sh

MODEL=M20
SETB=set1
SETE=set26

NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`
if [ $MODEL == "M00"  ];then
	rm -rf optimize/MODEL_M00; mkdir -p optimize/MODEL_M00
	\cp -r model_initial/* optimize/MODEL_M00
fi

# change parfile
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

mpirun -np $NPROC $fksem/bin/xfwat2_postproc_opt $MODEL $SETB $SETE true

for step in `grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`;do
  echo "======================================="
  echo  Meshing for model $MODEL step:$step
  echo "========================================"
  #=====
  LOCAL_PATH="./optimize/MODEL_${MODEL}_step${step}"
  ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH DATA/Par_file
  ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH  DATA/meshfem3D_files/Mesh_Par_file
  ./utils/change_par_file.sh SAVE_MESH_FILES '.false.' DATA/Par_file
  #sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${MODEL}_step${step}" DATA/Par_file
  #sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${MODEL}_step${step}" DATA/meshfem3D_files/Mesh_Par_file
  #sed -i '/SAVE_MESH_FILES                 =/c\SAVE_MESH_FILES                 = .false.' DATA/Par_file
 
  mpirun -np $NPROC $fksem/bin/xmeshfem3D
  mpirun -np $NPROC $fksem/bin/xgenerate_databases
done



