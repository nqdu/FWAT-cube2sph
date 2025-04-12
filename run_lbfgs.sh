#!/bin/bash
set_fwat1()
{
  local simu_type=$1
  local njobs=$2
  local start_set=$3
  local nevts=`cat src_rec/sources.dat.$simu_type|wc -l`
  local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`

  if [[ $iter == 0 || $iter == $iter_end ]]; then
    sed -i "/SHOW_DETAILS:/c\SHOW_DETAILS: .true." fwat_params/FWAT.PAR
  fi

  # substitute 
  cp sbash_fwat1_measure.sh tmp.fwat1.$simu_type.sh
  \cp DATA/Par_file.$simu_type DATA/Par_file
  local fwd=tmp.fwat1.$simu_type.sh
  sed -i "/#SBATCH --array=/c\#SBATCH --array=1-$narray%5" $fwd
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/NJOBS=/c\NJOBS=$njobs" $fwd
  sed -i "/START_SET=/c\START_SET=$start_set" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd
  sed -i "/DO_LS=/c\DO_LS=0" $fwd

  if [[ $simu_type == "noise" ]]; then
    sed -i "/#SBATCH --time=/c\#SBATCH --time=01:30:00" $fwd
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
  local nevts=`cat src_rec/sources.dat.ls.$simu_type|wc -l`
  local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`

  # substitute 
  cp sbash_fwat1_measure.sh tmp.fwat3.$simu_type.sh
  \cp DATA/Par_file.$simu_type DATA/Par_file
  local fwd=tmp.fwat3.$simu_type.sh
  sed -i "/#SBATCH --array=/c\#SBATCH --array=1-$narray%5" $fwd
  sed -i "/#SBATCH --job-name=/c\#SBATCH --job-name=LS" $fwd 
  sed -i "/#SBATCH --output=/c\#SBATCH --output=LS-%j_set%a.txt" $fwd 
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/NJOBS=/c\NJOBS=$njobs" $fwd
  sed -i "/DO_LS=/c\DO_LS=1" $fwd
  sed -i "/START_SET=/c\START_SET=$start_set" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd

  if [[ $simu_type == "noise" ]]; then
    sed -i "/#SBATCH --time=/c\#SBATCH --time=01:30:00" $fwd
  else
    sed -i "/#SBATCH --time=/c\#SBATCH --time=00:25:00" $fwd
  fi

  # run LS
  echo "LS $simu_type simulation for updated model ..."
}

# simu_type
simu_type=tele
NJOBS=1

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
sete=`cat src_rec/sources.dat.$simu_type |wc -l`
sete=`echo "$setb + $sete - 1" |bc`

# current model
iter=`grep ^iter: lbfgs.yaml | cut -d: -f2`
flag=`grep ^flag: lbfgs.yaml | cut -d: -f2`
echo $iter 
mod=M`printf %02d $iter`
echo "iteration $iter $mod $flag"

# check flag type and run 
if [ $flag == "INIT" ]; then 
  set_fwat1 $simu_type $NJOBS $setb
  job_adj=$(sbatch tmp.fwat1.$simu_type.sh|cut -d ' ' -f4 )
  
  # get lbfgs-direction
  fwd=sbash_fwat2_lbfgs.sh
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd
  sbatch --dependency=afterok:${job_adj} $fwd  
  
elif [ $flag == "GRAD"  ];then 

else  # line search
  #echo "set_fwat3 $simu_type $NJOBS $setb"
  #exit 
  set_fwat3 $simu_type $NJOBS $setb 
  job_line=$(sbatch tmp.fwat3.$simu_type.sh|cut -d ' ' -f4 )

  # 
  # check wolfe condition
  fwd=sbash_fwat3_wolfe.sh
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd
  sed -i "/MODEL=/c\MODEL=$mod" $fwd
  sbatch --dependency=afterok:${job_line} $fwd 
fi