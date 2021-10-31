#!/bin/bash
#SBATCH --nodes=4 
#SBATCH --ntasks-per-node=20
#SBATCH --time=01:59:59
#SBATCH --job-name FWD22
#SBATCH --output=FWD22_%j.txt
#SBATCH --partition=TH_HPC3
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
#SBATCH --mail-type=ALL

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
#module load gcc/9.2.0
#module load intel openmpi
#module load NiaEnv/2018a
#module load intel openmpi
#cd $SLURM_SUBMIT_DIR

#==== Comment out the following if running SEM mesh with new models====#
#mpirun $SEM3DROOT/bin/xmeshfem3D
#mpirun $SEM3DROOT/bin/xgenerate_databases
#exit
#======================================================================#

MODEL=M00
SET=set22
set_num=`echo $SET |awk -Fset '{printf"%d\n",$2}'`
fksem='/THL8/home/iggluyf/nqdu/specfem3d-joint/'
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`

# set_num > 10 for tele
if [ $set_num -le 10 ];then
  #mpirun -oversubscribe -np $NPROC $fksem/bin/xfwat1_fwd_measure_adj $MODEL $SET noise 3
  mpirun -np $NPROC $fksem/bin/xfwat0_forward_data $MODEL $SET noise
else
  #mpirun -oversubscribe -np $NPROC $fksem/bin/xfwat1_fwd_measure_adj $MODEL $SET tele 3
  mpirun  -np $NPROC $fksem/bin/xfwat0_forward_data $MODEL $SET tele
fi  

