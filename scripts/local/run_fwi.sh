#!/bin/bash
set -e

SET_FWAT()
{
  # local vars
  local flag=$1
  local mod=$2

  nsimtypes="${#SIMU_TYPES[@]}"
  for ((isim=0;isim<$nsimtypes;isim++)); 
  do 
    local simu_type=${SIMU_TYPES[$isim]}
    local nevts=`awk 'END { print NR }' src_rec/sources.dat.$simu_type`
    local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`
    local njobs=$nevts

    # copy files
    fwd=tmp.adj.$simu_type.sh
    \cp sbash_measure.sh $fwd

    # substitute 
    \cp DATA/Par_file.$simu_type DATA/Par_file
    sed -i "/NJOBS=/c\NJOBS=$njobs" $fwd
    sed -i "/simu_type=/c\simu_type=$simu_type" $fwd
    sed -i "/LOCAL_PC=/c\LOCAL_PC=1" $fwd

    # run forward/adjoint simulation
    echo "forward/adjoint $simu_type simulation for new model  ..."

    #  run simulation
    if [ "$flag" == "INIT" ]; then 
      bash $fwd > LOG/FWD_ADJ.$simu_type.$iter.txt
    else
      bash $fwd > LOG/LS.$simu_type.$iter.txt
    fi
  done 
}

set -e 
. parameters.sh

######### USER PARAMETERS
# parameters
####### STOP HERE #########

# mkdir 
mkdir -p misfits optimize solver LOG

# some jobid 
job_adj=0
job_post=0
job_step=0
job_line=0
job_wait=0

date
for ii in `seq 1 4`;do 

  # current model
  iter=`fwat-utils getparam iter fwat_params/lbfgs.yaml`
  flag=`fwat-utils getparam flag fwat_params/lbfgs.yaml`
  mod=M`printf %02d $iter`
  echo "iteration $iter $mod $flag"


  # copy first model to MODEL_M00
  if [[ "$iter" -eq 0 && ! -d "optimize/MODEL_M00" ]]; then
  #if [[ "$iter" -eq 0 ]]; then
    # set model name
    mkdir -p optimize/MODEL_M00
    \cp ./DATABASES_MPI/* optimize/MODEL_M00
  fi

  # check flag type and run 
  if [ $flag == "INIT" ]; then 
    SET_FWAT $flag $mod
    
    # sum kernels, get search direction, generate trial model 
    fwd=sbash_postproc_kl.sh
    bash $fwd > LOG/POST.$iter.txt 
    
  elif [ $flag == "GRAD"  ];then 
    # get search direction, generate trial model 
    fwd=sbash_postproc_kl.sh
    echo "Post processing ..."
    bash $fwd > LOG/POST.$iter.txt 

  else  # line search
    SET_FWAT $flag $mod

    # check wolfe condition
    fwd=sbash_wolfe.sh
    echo "checking wolfe condition ..."
    bash $fwd > LOG/WOLFE.$iter.txt 
  fi
done

date 