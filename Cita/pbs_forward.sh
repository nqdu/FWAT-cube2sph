#!/bin/bash -l
#PBS -l nodes=2:ppn=80
#PBS -l walltime=01:00:00
#PBS -N FWD22
#PBS -q starq

cd $PBS_O_WORKDIR
module load mpi/gcc-openmpi blas

#==== Comment out the following if running SEM mesh with new models====#
#mpirun $SEM3DROOT/bin/xmeshfem3D
#mpirun $SEM3DROOT/bin/xgenerate_databases
#exit
#======================================================================#

MODEL=M00
SET=set28
simu_type=rf
set_num=`echo $SET |awk -Fset '{printf"%d\n",$2}'`
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`

. utils/parameters.sh
mpirun -np $NPROC $fksem/bin/xfwat0_forward_data $MODEL $SET $simu_type



