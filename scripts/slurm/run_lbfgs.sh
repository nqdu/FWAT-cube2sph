#!/bin/bash
set -e

SET_FWAT()
{
  # local vars
  local njobs=$1
  local flag=$2
  local mod=$3

  local job_ids=()

  nsimtypes="${#SIMU_TYPES[@]}"
  for ((isim=0;isim<$nsimtypes;isim++)); 
  do 
    local simu_type=${SIMU_TYPES[$isim]}
    local nevts=`awk 'END { print NR }' src_rec/sources.dat.$simu_type`
    local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`

    # copy files
    fwd=tmp.adj.$simu_type.sh
    \cp sbash_measure.sh $fwd

    # substitute 
    \cp DATA/Par_file.$simu_type DATA/Par_file
    sed -i "/#SBATCH --array=/c\#SBATCH --array=1-$narray%5" $fwd
    sed -i "/MODEL=/c\MODEL=${mod}" $fwd
    sed -i "/NJOBS=/c\NJOBS=$njobs" $fwd
    sed -i "/START_SET=/c\START_SET=1" $fwd
    sed -i "/simu_type=/c\simu_type=$simu_type" $fwd
    sed -i "/LOCAL_PC=/c\LOCAL_PC=0" $fwd
    if [ "$flag" == "INIT" ]; then 
      sed -i "/#SBATCH --job-name=/c\#SBATCH --job-name=FWD_ADJ" $fwd
      #sed -i "/#SBATCH --output=/c\#SBATCH --output=FWD_ADJ-%j_set%a.txt" $fwd
      sed -i "/#SBATCH --output=/c\#SBATCH --output=LOG/FWD_ADJ-%j_set%a.txt" $fwd
    else
      sed -i "/#SBATCH --job-name=/c\#SBATCH --job-name=LS" $fwd
      #sed -i "/#SBATCH --output=/c\#SBATCH --output=LS-%j_set%a.txt" $fwd
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

######### USER PARAMETERS
# parameters
NJOBS=1
####### STOP HERE #########

source parameters.sh

# mkdir 
mkdir -p misfits optimize solver LOG

# some jobid 
job_adj=0
job_post=0

for ii in `seq 1 4`;do 

  # current model
  iter=`python $FWATLIB/get_param.py iter $FWATPARAM/lbfgs.yaml`
  flag=`python $FWATLIB/get_param.py flag $FWATPARAM/lbfgs.yaml`
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
    SET_FWAT $NJOBS $flag $mod
    
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
    SET_FWAT $NJOBS $flag $mod

    # check wolfe condition
    fwd=sbash_wolfe.sh
    #sed -i "/#SBATCH --output=/c\#SBATCH --output=LOG//WOLFE_%j.txt" $fwd
    echo "checking wolfe condition ..."
    job_post=$(sbatch --dependency=afterok:${job_adj} $fwd | cut -d ' ' -f4)
  fi

  # wait to finish
  srun --dependency=afterok:${job_post} --nodes=1 --time=00:01:05 --ntasks=1 --job-name=wait  ./wait.sh
done