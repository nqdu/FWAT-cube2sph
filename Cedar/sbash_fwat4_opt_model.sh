#!/bin/bash
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=24
#SBATCH --mem=0
#SBATCH --time=00:45:59
#SBATCH --job-name=OPT
#SBATCH --output=OPT_%j.txt

set -e 
#source activate pygmt 
module load intel openmpi

#####################################
###### Parameters ##################

# simu_type
simu_type=noise

# start model
MODEL=M20
SETB=set1
SETE=set26

# include file
. utils/parameters.sh

# go with tele data 
#startidx=$(printf %d `echo $startmod | cut -d'M' -f2`)
startidx=`echo $MODEL | cut -d'M' -f2`
jj=$(printf %02d `echo "$startidx + 0" | bc`)
echo M$jj

# compute misfits for each step size
cd plots
cd ..

# compute best model by using quadratic interpolation 
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`
mpirun -np $NPROC $fksem/bin/xfwat4_opt_model $MODEL $SETB $SETE true

# mesh new model
jnext=$(printf %02d $(echo "1+$startidx"|bc))
LOCAL_PATH="./optimize/MODEL_M${jnext}"
./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH DATA/Par_file
./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH  DATA/meshfem3D_files/Mesh_Par_file
./utils/change_par_file.sh SAVE_MESH_FILES '.false.' DATA/Par_file
mpirun -np $NPROC $fksem/bin/xmeshfem3D
mpirun -np $NPROC $fksem/bin/xgenerate_databases

# find the min misfit
#jnext=$(printf %02d $(echo "1+$startidx"|bc))
misfit_cur=`head -1 plots/M${jj}.mis.avg | awk '{printf("%g",$2)}'`
misfit_next=`tail -1 plots/M${jj}.mis.avg | awk '{printf("%g",$2)}'`
echo "misfit for iteration $jj and $jnext: $misfit_cur $misfit_next" >> misfit.log 
