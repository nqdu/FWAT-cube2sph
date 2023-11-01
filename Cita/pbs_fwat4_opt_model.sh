#!/bin/bash -l
#PBS -l nodes=2:ppn=80
#PBS -l walltime=01:30:00
#PBS -N OPT
#PBS -q starq
#PBS -j oe

cd $PBS_O_WORKDIR

set -e 
source activate pygmt 
module load mpi/gcc-openmpi blas 

#####################################
###### Parameters ##################

# simu_type
simu_type=rf

# start model
MODEL=M10
SETB=set23
SETE=set28

# include file
. utils/parameters.sh

# go with tele data 
#startidx=$(printf %d `echo $startmod | cut -d'M' -f2`)
startidx=`echo $MODEL | cut -d'M' -f2`
jj=$(printf %02d `echo "$startidx + 0" | bc`)
echo M$jj

# compute misfits for each step size
cd plots
bash plot_misfit/plt_line_search.multiband.${simu_type}.bash M$jj
#if [ $simu_type == "tele"  ];then
#	bash plot_misfit/plt_line_search.multiband.tele.bash M$jj
#else 
#	bash plot_misfit/plt_line_search.multiband.ANAT.bash M$jj
#fi
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

# save log file
#mkdir -p LOGS/MODEL_$jj
#mv  *.o[0-9]* *.txt LOGS/MODEL_$jj -f
