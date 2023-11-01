#!/bin/bash
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:47:59
#SBATCH --job-name FWD26
#SBATCH --output=FWD26_%j.txt
#SBATCH --mem=0
#SBATCH --account=rrg-liuqy
##SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
##SBATCH --mail-type=ALL

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
SET=set26
simu_type=rf
set_num=`echo $SET |awk -Fset '{printf"%d\n",$2}'`
fksem="/home/l/liuqy/nqdu/scratch/specfem3d//"
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

