#!/bin/bash
set_fwat1()
{
  local simu_type=$1
  local njobs=$2
  local start_set=$3
  local nevts=`cat src_rec/sources.dat.$simu_type|wc -l`
  local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`

  # substitute 
  cp sbash_measure.sh tmp.fwat1.$simu_type.sh
  \cp DATA/Par_file.$simu_type DATA/Par_file
  local fwd=tmp.fwat1.$simu_type.sh
  sed -i "/#SBATCH --array=/c\#SBATCH --array=1-$narray%5" $fwd
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/NJOBS=/c\NJOBS=$njobs" $fwd
  sed -i "/START_SET=/c\START_SET=$start_set" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd
  sed -i "/LOCAL_PC=/c\LOCAL_PC=0" $fwd

  if [[ $simu_type == "noise" ]]; then
    sed -i "/#SBATCH --time=/c\#SBATCH --time=00:25:00" $fwd
  else
    sed -i "/#SBATCH --time=/c\#SBATCH --time=00:25:00" $fwd
  fi

  # run forward/adjoint simulation
  echo "forward/adjoint $simu_type simulation for new model  ..."
}

#!/bin/bash
set_fwat3()
{
  local simu_type=$1
  local njobs=$2
  local start_set=$3
  local nevts=`cat src_rec/sources.dat.$simu_type|wc -l`
  local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`

  # substitute 
  cp sbash_measure.sh tmp.fwat3.$simu_type.sh
  \cp DATA/Par_file.$simu_type DATA/Par_file
  local fwd=tmp.fwat3.$simu_type.sh
  sed -i "/#SBATCH --array=/c\#SBATCH --array=1-$narray%5" $fwd
  sed -i "/#SBATCH --job-name=/c\#SBATCH --job-name=LS" $fwd
  sed -i "/#SBATCH --output=/c\#SBATCH --output=LS-%j_set%a.txt" $fwd
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/NJOBS=/c\NJOBS=$njobs" $fwd
  sed -i "/START_SET=/c\START_SET=$start_set" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd
  sed -i "/LOCAL_PC=/c\LOCAL_PC=0" $fwd


  if [[ $simu_type == "noise" ]]; then
    sed -i "/#SBATCH --time=/c\#SBATCH --time=00:25:00" $fwd
  else
    sed -i "/#SBATCH --time=/c\#SBATCH --time=00:25:00" $fwd
  fi

  # run forward/adjoint simulation
  echo "forward/adjoint $simu_type simulation for new model  ..."
}

set -e 
. parameters.sh

# simu_type
simu_type=noise
NJOBS=8

# mkdir 
mkdir -p misfits optimize solver

# some jobid 
job_adj=0
job_post=0
job_step=0
job_line=0
job_wait=0

# set number
# initialize set number
setb=1

for ii in `seq 1 4`;do 

  # current model
  iter=`python $FWATLIB/get_param.py iter $FWATPARAM/lbfgs.yaml`
  flag=`python $FWATLIB/get_param.py flag $FWATPARAM/lbfgs.yaml`
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
    set_fwat1 $simu_type $NJOBS $setb
    job_adj=$(sbatch tmp.fwat1.$simu_type.sh|cut -d ' ' -f4 )
    
    # sum kernels, get search direction, generate trial model 
    fwd=sbash_postproc_kl.sh
    sed -i "/SOURCE_FILE=/c\SOURCE_FILE=./src_rec/sources.dat.$simu_type" $fwd
    job_post=$(sbatch --dependency=afterok:${job_adj} $fwd | cut -d ' ' -f4)
    
  elif [ $flag == "GRAD"  ];then 
    # get search direction, generate trial model 
    fwd=sbash_postproc_kl.sh
    sed -i "/SOURCE_FILE=/c\SOURCE_FILE=./src_rec/sources.dat.$simu_type" $fwd
    echo "Post processing ..."
    job_post=$(sbatch $fwd | cut -d ' ' -f4)

  else  # line search
    set_fwat3 $simu_type $NJOBS $setb 
    job_line=$(sbatch tmp.fwat3.$simu_type.sh|cut -d ' ' -f4 )

    # check wolfe condition
    fwd=sbash_wolfe.sh
    sed -i "/SIMU_TYPE=/c\SIMU_TYPE=$simu_type" $fwd
    sed -i "/SOURCE_FILE=/c\SOURCE_FILE=./src_rec/sources.dat.$simu_type" $fwd
    echo "checking wolfe condition ..."
    job_post=$(sbatch --dependency=afterok:${job_line} $fwd | cut -d ' ' -f4)
  fi

  # wait to finish
  srun --dependency=afterok:${job_post} --nodes=1 --time=00:00:10 --ntasks=1 --job-name=wait  ./wait.sh
done