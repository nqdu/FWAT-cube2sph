#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=160
#SBATCH --time=00:16:59
#SBATCH --job-name OPT
#SBATCH --output=OPT_%j.txt
#SBATCH --partition=compute

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
set -e 

. parameters.sh 
source module_env

NPROC=$SLURM_NTASKS
model=M24
SIMU_TYPE=tele

# get step_fac/dmax
MIS_FILE=misfits/$model.mis
step_fac=`tail -1 $MIS_FILE |cut -d'=' -f2 |awk '{print $1}'`
dmax=`tail -1 $MIS_FILE |cut -d'=' -f2 |awk '{print $2}'`

# calculate misfit
#cd plots
info=`python cal_misfit.py $model $SIMU_TYPE 00`
chi=`echo $info |awk '{print $1}'`
info=`python cal_misfit.py $model $SIMU_TYPE 01`
chi1=`echo $info |awk '{print $1}'`
#cd ..
echo $chi1 $chi 
echo " " >> $MIS_FILE
echo "misfit = $chi $chi1" >> $MIS_FILE
#echo "misfit for model $model and next  = $chi $chi1" >> misfit.log

# check if wolfe condition is satisfied 
echo "python $OPT_LIB/wolfe_cond.py $model $SIMU_TYPE $chi $chi1 $step_fac $NPROC"
echo " " >> $MIS_FILE
python $OPT_LIB/wolfe_cond.py $model $SIMU_TYPE $chi $chi1 $step_fac $NPROC >> $MIS_FILE
step_fac_opt=`tail -1 $MIS_FILE |cut -d'=' -f2 |awk '{print $1}'`
echo $step_fac_opt
LSDIR=./optimize/MODEL_${model}_step01
if [[ `grep "failed" $MIS_FILE |wc -l` == 1 ]]; then  
  echo "resubmit this job by using step_fac = $step_fac_opt" 
  python $FWATLIB/set_param.py STEP_FAC $step_fac_opt fwat_params/lbfgs.yaml 
  exit 1
else 
  chimin=`tail -1 $MIS_FILE |cut -d'=' -f2 |awk '{print $2}'`
  if (( `echo "$step_fac $step_fac_opt" |awk '{print sqrt(($1-$2)^2)  >= 0.001 }'`  )); then 
    step_fac=$step_fac_opt
    python $OPT_LIB/get_lbfgs_step_fac.py $model $LSDIR $step_fac $NPROC >> $MIS_FILE

    # regenerate database
    LOCAL_PATH=./DATABASES_MPI
    $change_par LOCAL_PATH $LSDIR DATA/Par_file
    $change_par LOCAL_PATH $LSDIR  DATA/meshfem3D_files/Mesh_Par_file
    $change_par SAVE_MESH_FILES .false. DATA/Par_file
    echo -e ".false.\n.true." > adepml_stage
    mpirun -np $NPROC $fksem/bin/xgenerate_databases
    \rm adepml_stage
  fi
  
  per=`echo $step_fac $dmax |awk '{printf "%g", $1*$2}'`
  echo "optimial step_fac = $step_fac percent = $per"
  icur=$(echo $model |awk -F'M' '{print $2}')
  inext=$(printf "%02d" `echo $model |awk -F'M' '{print $2+1}'`)
  echo misfit for iteration $icur and $inext $chi $chimin >> misfit.log
  
  rm -rf ./optimize/MODEL_M$inext 
  mv $LSDIR ./optimize/MODEL_M$inext 
fi 

fwd=output_fwat4_log_${model}.txt
echo "******************************************************" > $fwd
echo " Finished FWAT stage4 here!!!" >> $fwd
echo " " >> $fwd
