#!/bin/bash -l
#PBS -l nodes=2:ppn=80
#PBS -l walltime=01:30:00
#PBS -N FWD_ADJ
#PBS -J 23-28
#PBS -q starq
#PBS -j oe

cd $PBS_O_WORKDIR

module load mpi/gcc-openmpi  blas

# include file
. utils/parameters.sh

#==== Comment out the following if running SEM mesh with new models====#
MODEL=M10
simu_type=rf

# SET
SET=set$PBS_ARRAY_INDEX
set_num=`echo $SET |awk -Fset '{printf"%d\n",$2}'`
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`

# set_num > 10 for tele
mpirun -np $NPROC $fksem/bin/xfwat1_fwd_measure_adj $MODEL $SET $simu_type 3
