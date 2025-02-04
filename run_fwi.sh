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
  cp sbash_fwat1_fwd_measure_adj.sh tmp.fwat1.$simu_type.sh
  \cp DATA/Par_file.$simu_type DATA/Par_file
  local fwd=tmp.fwat1.$simu_type.sh
  sed -i "/#SBATCH --array=/c\#SBATCH --array=1-$narray%5" $fwd
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/NJOBS=/c\NJOBS=$njobs" $fwd
  sed -i "/START_SET=/c\START_SET=$start_set" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd

  if [[ $simu_type == "noise" ]]; then
    sed -i "/#SBATCH --time=/c\#SBATCH --time=01:30:00" $fwd
  else
    sed -i "/#SBATCH --time=/c\#SBATCH --time=00:35:00" $fwd
  fi

  # run forward/adjoint simulation
  echo "computing forward/adjoint simulation for $simu_type ..."
}

#!/bin/bash
set_fwat3()
{
  local simu_type=$1
  local njobs=$2
  local start_set=$3
  local nevts=`cat src_rec/sources.dat.$simu_type|wc -l`

  # njobs *2
  njobs=`echo "$njobs *2"|bc`
  local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`

  # substitute 
  cp sbash_fwat3_linesearch.sh tmp.fwat3.$simu_type.sh
  \cp DATA/Par_file.$simu_type DATA/Par_file
  local fwd=tmp.fwat3.$simu_type.sh
  sed -i "/#SBATCH --array=/c\#SBATCH --array=1-$narray%5" $fwd
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/NJOBS=/c\NJOBS=$njobs" $fwd
  sed -i "/START_SET=/c\START_SET=$start_set" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd

  if [[ $simu_type == "noise" ]]; then
    sed -i "/#SBATCH --time=/c\#SBATCH --time=01:30:00" $fwd
  else
    sed -i "/#SBATCH --time=/c\#SBATCH --time=00:35:00" $fwd
  fi

  # run forward/adjoint simulation
  echo "computing line search for $simu_type ..."
}

# parameters
iter_start=55
iter_end=55

# L-BFGS params
lbfgs_start=38

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

for iter in `seq $iter_start $iter_end`;do
  # current model
  mod=M`printf %02d $iter`
  mod_lbfgs_start=M`printf %02d $lbfgs_start`
  echo "iteration $iter $mod $mod_lbfgs_start"

  # create misfit file 
  :> plots/$mod.mis

  # run RF/tele simulation
  set_fwat1 $simu_type $NJOBS $setb
  #exit
  job_adj=$(sbatch tmp.fwat1.$simu_type.sh|cut -d ' ' -f4 )

  # run line search
  if [ ${mod} == $mod_lbfgs_start  ]; then 
    echo ${mod} > lbfgs.in
    echo "-1" >> lbfgs.in
  else 
    info=`head -1 lbfgs.in`
    echo $info > lbfgs.in 
    echo "1" >> lbfgs.in
  fi
  #exit

  # run postprocessing
  fwd=sbash_fwat2_postproc_opt.sh 
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/SOURCE_FILE=/c\SOURCE_FILE=./src_rec/sources.dat.$simu_type" $fwd
  sed -i "/SETE=/c\SETE=${sete}" $fwd
  \cp fwat_params/FWAT.PAR.$simu_type  fwat_params/FWAT.PAR
  echo "postprocessing ..."
  
  job_post=$(sbatch --dependency=afterok:${job_adj} $fwd | cut -d ' ' -f4)

  # run line search
  set_fwat3 $simu_type $NJOBS $setb
  job_line=$(sbatch --dependency=afterok:${job_post} tmp.fwat3.$simu_type.sh | cut -d ' ' -f4)

  # generate next model
  echo "generating opt model ..."
  fwd=sbash_fwat4_opt_model.sh 
  sed -i "/model=/c\model=${mod}" $fwd
  sed -i "/SIMU_TYPE=/c\SIMU_TYPE=${simu_type}" $fwd
  job_step=$(sbatch --dependency=afterok:${job_line} $fwd | cut -d ' ' -f4)

  # wait job to finish
  fwd=wait.sh
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  srun --dependency=afterok:${job_step} --partition=compute --nodes=1 --ntasks=1 --time=00:15:02 wait.sh 
  mkdir -p LOG/$mod
  mv *.txt LOG/$mod
done
