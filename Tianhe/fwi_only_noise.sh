#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=20
#SBATCH --job-name FWI
#SBATCH --output=FWI_%j.log
#SBATCH --partition=TH_HPC3

set -e 

#####################################
###### Parameters ##################

# event sets
setb=1
sete=10

# simu_type
simu_type=noise 

# start model
startmod=M00

# iterations 
iters=5

# HPC facilities
nnodes=8
ntasks_per_node=20

# FKSEM path
fksem='/THL8/home/iggluyf/nqdu/specfem3d-joint/'

###### Parameters END ################
####################################

# LOAD all module
source activate base

echo "========================="
echo "FWI begin `date`"
echo "========================="

# copy parameters if required
if [ $simu_type == 'noise' ];then
    echo "copying params for noise"
    cp DATA/Par_file.noise DATA/Par_file
    cp fwat_params/FWAT.PAR.noise fwat_params/FWAT.PAR
    cp fwat_params/MEASUREMENT.PAR.noise fwat_params/MEASUREMENT.PAR
    cp src_rec/sources_ls.dat.noise src_rec/sources_ls.dat
elif [ $simu_type == 'tele' ];then
    echo "copying params for tele"
    cp DATA/Par_file.tele DATA/Par_file
    cp fwat_params/FWAT.PAR.tele fwat_params/FWAT.PAR
    cp fwat_params/MEASUREMENT.PAR.tele fwat_params/MEASUREMENT.PAR
    cp src_rec/sources_ls.dat.tele src_rec/sources_ls.dat
fi 

# machinefile 
NPROC=`grep "^NPROC" DATA/Par_file | awk '{print $3}'`
srun hostname | sort -n | uniq | awk -F. '{print $1}' > slurm.host 
python utils/generate_hostfile.py $nnodes $ntasks_per_node $NPROC
nodes=`ls slurm.host.? | wc -l`

# generate mesh if required
if [ $startmod == "M00" ]; then 
    LOCAL_PATH="OUTPUT_FILES/DATABASES_MPI/"
    rm -rf $LOCAL_PATH/*.bin

    # initial model
    LOCAL_PATH="OUTPUT_FILES/DATABASES_MPI/"
    mkdir -p $LOCAL_PATH
    rm -rf $LOCAL_PATH/*.bin
    ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH DATA/Par_file
    ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH  DATA/meshfem3D_files/Mesh_Par_file
    cp initial_model/* ${LOCAL_PATH}/ -r 
    echo 'mpirun --oversubscribe -np $NPROC $fksem/bin/xmeshfem3D'
    mpirun -machinefile slurm.host -np $NPROC $fksem/bin/xmeshfem3D
    echo 'mpirun --oversubscribe -np $NPROC $fksem/bin/xgenerate_databases'
    mpirun -machinefile slurm.host -np $NPROC $fksem/bin/xgenerate_databases

fi

# go with tele data 
startidx=$(printf %d `echo $startmod | cut -d'M' -f2`)
for((ii=0;ii<$iters;ii++));
do  
    echo " "
    echo "==============================="
    echo "Iteration $ii"
    echo "==============================="
    
    # current and next index
    jj=$(printf %02d $(echo "$ii+$startidx"|bc))
    jnext=$(printf %02d $(echo "$ii+1+$startidx"|bc))

    # change LOCAL_PATH to the model database
    if [ $jj == "00" ];then
        ./utils/change_par_file.sh LOCAL_PATH ./OUTPUT_FILES/DATABASES_MPI DATA/Par_file
        ./utils/change_par_file.sh LOCAL_PATH ./OUTPUT_FILES/DATABASES_MPI DATA/meshfem3D_files/Mesh_Par_file
    else
        ./utils/change_par_file.sh LOCAL_PATH ./optimize/MODEL_M${jj} DATA/Par_file
        ./utils/change_par_file.sh LOCAL_PATH ./optimize/MODEL_M${jj} DATA/meshfem3D_files/Mesh_Par_file
    fi 
    
    # forward and adjoint simulation
    echo "Forward and Adjoint simulations ..."
    python src/run_fwdadj.py M$jj $setb $sete $simu_type $NPROC $nodes 

    # # post-processing
    echo "Post processing ..."
    mkdir -p optimize
    if [ $jj == "00" ]; then 
        rm -rf optimize/MODEL_M00
        cp -r initial_model optimize/MODEL_M00
    fi 
    echo "mpirun -machinefile slurm.host -np $NPROC $fksem/bin/xfwat2_postproc_opt M$jj set$setb set$sete true"
    mpirun -machinefile slurm.host -np $NPROC $fksem/bin/xfwat2_postproc_opt M$jj set$setb set$sete true > POST.txt 
    for step in `grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`;do
        echo "======================================="
        echo  Meshing for model $MODEL step:$step
        echo "========================================"
        #=====
        ./utils/change_par_file.sh LOCAL_PATH ./optimize/MODEL_M${jj}_step${step} DATA/Par_file
        ./utils/change_par_file.sh LOCAL_PATH ./optimize/MODEL_M${jj}_step${step} DATA/meshfem3D_files/Mesh_Par_file
        ./utils/change_par_file.sh SAVE_MESH_FILES .false. DATA/Par_file       
        mpirun -machinefile slurm.host -np $NPROC $fksem/bin/xmeshfem3D
        mpirun -machinefile slurm.host -np $NPROC $fksem/bin/xgenerate_databases
    done

    # linesearch
    echo "Line Searching ..."
    if [[ $simu_type == "tele" ]];then 
        cp src_rec/sources_ls.dat.tele src_rec/sources_ls -r 
    else
        cp src_rec/sources_ls.dat.noise src_rec/sources_ls -r 
    fi 
    python src/run_linesearch.py M$jj $simu_type $NPROC $nodes
    
    # compute misfits for each step size
    cd plots 
    if [ $simu_type == "tele" ]; then 
        bash plot_misfit/plt_line_search.multiband.tele.bash M$jj
    else 
        bash plot_misfit/plt_line_search.multiband.ANAT.bash M$jj
    fi 

    # find the min misfit and copy to Model_next
    minval=`tail -n +2 M${jj}.mis.avg| awk '{print $2}' | sort -n | head -1`
    step=`grep $minval M${jj}.mis.avg -n | cut -d: -f2 | awk '{print $1}'`
    cd ../optimize
    cp MODEL_M${jj}_step${step} MODEL_M$jnext -rf 
    cd ..
    echo "misfit for iteration $ii with step_size $step : $minval"

    # save log file
    mkdir -p LOGS/MODEL_$jj
    mv  *.txt LOGS/MODEL_$jj -f

    echo " "

done 

echo "========================="
echo "FWI finished `date`"
echo "========================="