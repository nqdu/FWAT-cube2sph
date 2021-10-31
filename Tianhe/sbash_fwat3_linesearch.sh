#!/bin/bash
#SBATCH --nodes=4 
#SBATCH --ntasks-per-node=20
#SBATCH --time=01:59:59
#SBATCH --job-name LS0.050
#SBATCH --output=LS0.050_%j.txt
#SBATCH --partition=TH_HPC3
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
#SBATCH --mail-type=ALL

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
#module load intel/15.0.2 openmpi/intel/1.6.4
#module load intel openmpi
#module load NiaEnv/2018a
#module load intel openmpi 

#=====
#cd $PBS_O_WORKDIR
#cd $SLURM_SUBMIT_DIR
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/intel/lib/intel64:/usr/local/openmpi/lib
#=====
fksem='/THL8/home/iggluyf/nqdu/specfem3d-joint/'
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`
MODEL=M05_step0.050
SET=ls
SIMU_TYPE=tele
mpirun -np $NPROC $fksem/bin/xfwat3_linesearch $MODEL $SET $SIMU_TYPE

