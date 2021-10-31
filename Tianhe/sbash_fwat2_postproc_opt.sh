#!/bin/bash
#SBATCH --nodes=4 
#SBATCH --ntasks-per-node=20
#SBATCH --time=01:59:59
#SBATCH --job-name=POST
#SBATCH --output=POST_%j.txt
#SBATCH --partition=TH_HPC3

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
#module load intel/15.0.2 openmpi/intel/1.6.4
#module load intel openmpi
#module load NiaEnv/2018a
#module load intel  openmpi

#=====
#cd $PBS_O_WORKDIR
#cd $SLURM_SUBMIT_DIR
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/intel/lib/intel64:/usr/local/openmpi/lib
#=====

set -e
MODEL=M06
SETB=set12
SETE=set22
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`
fksem='/THL8/home/iggluyf/nqdu/specfem3d-joint/'
mpirun -np $NPROC $fksem/bin/xfwat2_postproc_opt $MODEL $SETB $SETE true

for step in `grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`;do
  echo "======================================="
  echo  Meshing for model $MODEL step:$step
  echo "========================================"
  #=====
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${MODEL}_step${step}" DATA/Par_file
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${MODEL}_step${step}" DATA/meshfem3D_files/Mesh_Par_file
  sed -i '/SAVE_MESH_FILES                 =/c\SAVE_MESH_FILES                 = .false.' DATA/Par_file
 
  mpirun -np $NPROC $fksem/bin/xmeshfem3D
  mpirun -np $NPROC $fksem/bin/xgenerate_databases
done



