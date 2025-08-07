#!/bin/bash
set -e

SET_FWAT()
{
  # local vars
  local flag=$1
  local mod=$2

  local job_ids=()

  nsimtypes="${#SIMU_TYPES[@]}"
  for ((isim=0;isim<$nsimtypes;isim++)); 
  do 
    local njobs=${NJOBS_PER_JOBARRAY[$isim]}
    local simu_type=${SIMU_TYPES[$isim]}
    local nevts=`awk 'END { print NR }' src_rec/sources.dat.$simu_type`
    local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`

    # copy files
    fwd=tmp.adj.$simu_type.sh
    \cp sbash_measure.sh $fwd

    # substitute 
    \cp DATA/Par_file.$simu_type DATA/Par_file
    sed -i "/#SBATCH --array=/c\#SBATCH --array=1-$narray%5" $fwd
    sed -i "/NJOBS=/c\NJOBS=$njobs" $fwd
    sed -i "/simu_type=/c\simu_type=$simu_type" $fwd
    sed -i "/LOCAL_PC=/c\LOCAL_PC=0" $fwd
    if [ "$flag" == "INIT" ]; then 
      sed -i "/#SBATCH --job-name=/c\#SBATCH --job-name=FWD_ADJ" $fwd
      sed -i "/#SBATCH --output=/c\#SBATCH --output=LOG/FWD_ADJ-%j_set%a.txt" $fwd
    else
      sed -i "/#SBATCH --job-name=/c\#SBATCH --job-name=LS" $fwd
      sed -i "/#SBATCH --output=/c\#SBATCH --output=LOG/LS-%j_set%a.txt" $fwd
    fi

    if [[ $simu_type == "noise" ]]; then
      sed -i "/#SBATCH --time=/c\#SBATCH --time=00:25:00" $fwd
    else
      sed -i "/#SBATCH --time=/c\#SBATCH --time=00:25:00" $fwd
    fi

    # run forward/adjoint simulation
    echo "forward/adjoint $simu_type simulation for new model  ..."

    # submit and get job id
    local jid=$(sbatch $fwd |cut -d ' ' -f4 ) 
    job_ids+=($jid)
  done 

  # return dependency strings
  local depend_string=$(IFS=:; echo "${job_ids[*]}")
  job_adj=$depend_string
}

######### USER PARAMETERS ###############
source parameters.sh

# mkdir 
mkdir -p misfits optimize solver LOG

# some jobid 
job_adj=0
job_post=0

for ii in `seq 1 4`;do 

  # current model
  iter=`fwat-utils getparam iter fwat_params/lbfgs.yaml`
  flag=`fwat-utils getparam flag fwat_params/lbfgs.yaml`
  mod=M`printf %02d $iter`
  #mkdir -p LOG/${mod}
  echo "iteration $iter $mod $flag"

  # copy first model to MODEL_M00
  if [[ "$iter" -eq 0 && ! -d "optimize/MODEL_M00" ]]; then
    mkdir -p optimize/MODEL_M00
    \cp ./DATABASES_MPI/* optimize/MODEL_M00
  fi

  # check flag type and run 
  if [ $flag == "INIT" ]; then 
    SET_FWAT $flag $mod
    
    # sum kernels, get search direction, generate trial model 
    fwd=sbash_postproc_kl.sh
    #sed -i "/#SBATCH --output=/c\#SBATCH --output=LOG/${mod}/POST_%j.txt" $fwd
    job_post=$(sbatch --dependency=afterok:${job_adj} $fwd | cut -d ' ' -f4)
    
  elif [ $flag == "GRAD"  ];then 
    # get search direction, generate trial model 
    fwd=sbash_postproc_kl.sh
    #sed -i "/#SBATCH --output=/c\#SBATCH --output=LOG/${mod}/POST_%j.txt" $fwd
    echo "Post processing ..."
    job_post=$(sbatch $fwd | cut -d ' ' -f4)

  else  # line search
    SET_FWAT $flag $mod

    # check wolfe condition
    fwd=sbash_wolfe.sh
    #sed -i "/#SBATCH --output=/c\#SBATCH --output=LOG//WOLFE_%j.txt" $fwd
    echo "checking wolfe condition ..."
    job_post=$(sbatch --dependency=afterok:${job_adj} $fwd | cut -d ' ' -f4)
  fi

  # wait to finish
  srun --dependency=afterok:${job_post} --nodes=1 --time=00:01:05 --ntasks=1 --job-name=wait  ./wait.sh
done