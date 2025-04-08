#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=160
#SBATCH --time=01:20:59
#SBATCH --job-name FWD22
#SBATCH --output=FWD22_%j.txt
#SBATCH --partition=compute

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
#module load gcc/9.2.0
#module load intel openmpi
module load NiaEnv/2019b
module load intel openmpi
#cd $SLURM_SUBMIT_DIR

#==== Comment out the following if running SEM mesh with new models====#
#mpirun $SEM3DROOT/bin/xmeshfem3D
#mpirun $SEM3DROOT/bin/xgenerate_databases
#exit
#======================================================================#

MODEL=M00
SET=set22
simu_type=noise
set_num=`echo $SET |awk -Fset '{printf"%d\n",$2}'`
fksem="/home/l/liuqy/nqdu/specfem3d//"
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`


mpirun -np $NPROC $fksem/bin/xfwat0_forward_data $MODEL $SET $simu_type
# set_num > 10 for tele
#if [ $set_num -le 22 ];then
#  #mpirun -oversubscribe -np $NPROC $fksem/bin/xfwat1_fwd_measure_adj $MODEL $SET noise 3
#  mpirun -oversubscribe -np $NPROC $fksem/bin/xfwat0_forward_data $MODEL $SET noise
#else
#  #mpirun -oversubscribe -np $NPROC $fksem/bin/xfwat1_fwd_measure_adj $MODEL $SET tele 3
#  mpirun -oversubscribe -np $NPROC $fksem/bin/xfwat0_forward_data $MODEL $SET tele
#fi  

