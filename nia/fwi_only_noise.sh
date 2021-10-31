#!/bin/bash
set -e 

#####################################
###### Parameters ##################

# event sets
setb=1
sete=10 

# start model
startmod=M00

# iterations 
iters=3

# partition
partition=debug 

###### Parameters END ################
####################################

# LOAD all module
source activate base
module load intel openmpi
. checkpoint.sh 

# change partition
sed -i "s/`grep '^partition=' submit_job_fwat1.sh`/partition=$partition/g" submit_job_fwat1.sh
sed -i "s/`grep '^partition=' submit_job_fwat2.sh`/partition=$partition/g" submit_job_fwat2.sh
sed -i "s/`grep '^partition=' submit_job_fwat3.sh`/partition=$partition/g" submit_job_fwat3.sh

echo "========================="
echo "FWI begin `date`"
echo "========================="

if [ $startmod == "M00" ]; then 
    LOCAL_PATH="OUTPUT_FILES/DATABASES_MPI/"
    rm -rf $LOCAL_PATH/*.bin

    # generate initial model
    fksem='/home/l/liuqy/nqdu/specfem3d/'
    NPROC=`grep "^NPROC" DATA/Par_file | awk '{print $3}'`
    cp initial_model/* ${LOCAL_PATH}/ -r 
    bash sbash_fwat0_mesh.sh 
    :>MESH
    # check finished
    fini=`ls | grep MESH |wc -l`
    while [ $fini -ne 1 ]
    do 
        sleep 5 
        fini=`ls | grep MESH |wc -l`
    done 
    sleep 30
    rm MESH*
fi 

# go with noise data 
startidx=$(printf %d `echo $startmod | cut -d'M' -f2`)
for((ii=$startidx;ii<$iters;ii++));
do  
    echo "==============================="
    echo "Iteration $ii"
    echo "==============================="
    # forward and adjoint simulation
    j=`printf %02d $ii`
    ./submit_job_fwat1.sh M$j noise $setb $sete 
    CheckFWDADJFinished M$j $setb $sete 

    # post-processing
    ./submit_job_fwat2.sh M$j $setb $sete  
    CheckPostFinished M$j 

    # linesearch
    ./submit_job_fwat3.sh M$j noise 
    CheckLineSearch M$j
    
    # compute misfits for each step size
    cd plots 
    bash plot_misfit/plt_line_search.multiband.ANAT.bash M$j 

    # next index
    jnext=$(printf %02d $(echo "$ii+1"|bc))

    # find the min misfit and copy to Model_next
    minval=`awk '{print $2}' M${j}.mis.avg | sort -n | head -1`
    step=`grep $minval M${j}.mis.avg -n | cut -d: -f2 | awk '{print $1}'`
    cd ../optimize
    cp MODEL_M${j}_step${step} MODEL_M$jnext -rf 
    cd ..
    echo "misfit for iteration $ii : $minval"

    # save log file
    mkdir -p LOGS/MODEL_$j
    mv output_fwat* FWD* POST*  LS* LOGS/MODEL_$j -f

    echo " "

done 

echo "========================="
echo "FWI finished `date`"
echo "========================="