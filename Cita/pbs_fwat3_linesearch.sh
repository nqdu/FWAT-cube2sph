#!/bin/bash -l
#PBS -l nodes=2:ppn=80
#PBS -l walltime=01:30:00
#PBS -N LS
#PBS -J 0-4
#PBS -q starq
#PBS -j oe 

cd $PBS_O_WORKDIR

module load mpi/gcc-openmpi blas 
. utils/parameters.sh 

NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`
sizes=(`grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`)
model=M10
MODEL=${model}_step${sizes[$PBS_ARRAY_INDEX]}
SET=ls
SIMU_TYPE=rf
mpirun -np $NPROC $fksem/bin/xfwat3_linesearch $MODEL $SET $SIMU_TYPE

