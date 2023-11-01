#!/bin/bash
set -e 

run_fwat1()
{
  local setb_tm=$1
  local sete_tm=$2
  local simu_type=$3 

  # copy files 
  if [ $simu_type == 'noise' ];then
    echo "copying params for noise"
    \cp DATA/Par_file.noise DATA/Par_file
    \cp fwat_params/FWAT.PAR.noise fwat_params/FWAT.PAR 
    \cp fwat_params/MEASUREMENT.PAR.noise fwat_params/MEASUREMENT.PAR
  elif [ $simu_type == 'tele' ];then
    echo "copying params for tele"
    \cp DATA/Par_file.tele DATA/Par_file
    \cp fwat_params/FWAT.PAR.tele fwat_params/FWAT.PAR
    \cp fwat_params/MEASUREMENT.PAR.tele fwat_params/MEASUREMENT.PAR
  elif [ $simu_type == "rf"  ];then
    echo "copying params for rf"
    \cp DATA/Par_file.rf DATA/Par_file
    \cp fwat_params/FWAT.PAR.rf fwat_params/FWAT.PAR
    \cp fwat_params/MEASUREMENT.PAR.rf fwat_params/MEASUREMENT.PAR
  else
    echo "$simu_type should be in [rf tele noise]"
    exit 1
  fi

  if [[ $iter == 0 || $iter == $iter_end ]]; then
    sed -i "/SHOW_DETAILS:/c\SHOW_DETAILS: .true." fwat_params/FWAT.PAR
  fi 

  # substitute 
  local fwd=sbash_fwat1_fwd_measure_adj.sh
  sed -i "/#SBATCH --array=/c\#SBATCH --array=$setb_tm-$sete_tm" $fwd
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd

  # run forward/adjoint simulation
  echo "computing forward/adjoint simulation ..."
  job_adj=$(sbatch $fwd |cut -d ' ' -f4 )
}

run_fwat3()
{
  local simu_type=$1

  if [ $simu_type == "tele" ];then
    \cp DATA/Par_file.tele DATA/Par_file
    \cp fwat_params/FWAT.PAR.tele fwat_params/FWAT.PAR
    \cp fwat_params/MEASUREMENT.PAR.tele fwat_params/MEASUREMENT.PAR
    \cp src_rec/sources_ls.dat.tele src_rec/sources_ls.dat -r 
  elif [ $simu_type == "noise"  ]; then
    \cp DATA/Par_file.noise DATA/Par_file
    \cp fwat_params/FWAT.PAR.noise fwat_params/FWAT.PAR 
    \cp fwat_params/MEASUREMENT.PAR.noise fwat_params/MEASUREMENT.PAR
    \cp src_rec/sources_ls.dat.noise src_rec/sources_ls.dat -r 
  else
    \cp DATA/Par_file.rf DATA/Par_file
    \cp fwat_params/FWAT.PAR.rf fwat_params/FWAT.PAR
    \cp fwat_params/MEASUREMENT.PAR.rf fwat_params/MEASUREMENT.PAR
    \cp src_rec/sources_ls.dat.rf src_rec/sources_ls.dat -r 
  fi

  local LOCAL_PATH="./OUTPUT_FILES/DATABASES_MPI"
  ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH DATA/Par_file
  ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH  DATA/meshfem3D_files/Mesh_Par_file

  local fwd=sbash_fwat3_linesearch.sh 
  local nsteps=`grep NUM_STEP fwat_params/FWAT.PAR |awk -F: '{print $2-1}'`
  sed -i "/#SBATCH --array=/c\#SBATCH --array=0-$nsteps" $fwd
  sed -i "/SIMU_TYPE=/c\SIMU_TYPE=$simu_type" $fwd
  sed -i "/model=/c\model=$mod" $fwd
  echo "line search ..."
  job_line=$(sbatch $fwd |cut -d ' ' -f4)
}

# parameters
iter_start=19
iter_end=20

# two dataset
simu_type1=noise
simu_type2=rf 
setb1=1
sete1=20
setb2=21
sete2=26

# make sure contiguous
setb=$setb1
sete=$sete2 

# mkdir 
mkdir -p misfits optimize solver 

# some jobid 
job_adj=0
job_post=0
job_line=0
job_misfit=0
job_wait=0

for iter in `seq $iter_start $iter_end`;do 
  # current model
  mod=M`printf %02d $iter`
  echo "iteration $iter $mod"

  # create misfit file 
  :> plots/$mod.mis.avg

  # run noise simulation
  run_fwat1 $setb1 $sete1 $simu_type1
  srun --dependency=afterok:${job_adj} --partition=compute --time=00:15:02 wait.sh

  # run RF/tele simulation
  run_fwat1 $setb2 $sete2 $simu_type2
  srun --dependency=afterok:${job_adj} --partition=compute --time=00:15:02 wait.sh 
  
  # save misfit for iter 0
  iter1=$iter
  iter=0
  if [[ "$iter" == "0" ]];then
    :> misfit0.log
    :> optimize/weight_kl.txt
    mis1=`bash step4_calmisfit.sh $iter $simu_type1 |cut -d' ' -f1`
    mis2=`bash step4_calmisfit.sh $iter $simu_type2 |cut -d' ' -f1` 
    for iset in `seq $setb1 $sete1`; do echo "1.0 / $mis1"|bc -l >> optimize/weight_kl.txt; done
    for iset in `seq $setb2 $sete2`; do echo "20.0 / $mis2"|bc -l >> optimize/weight_kl.txt; done 
    echo "$mis1 $mis2" > misfit0.log
  fi
  iter=$iter1

  # run postprocessing
  fwd=sbash_fwat2_postproc_opt.sh 
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/SETB=/c\SETB=set${setb}" $fwd
  sed -i "/SETE=/c\SETE=set${sete}" $fwd 
  echo "postprocessing ..."
  job_post=$(sbatch $fwd | cut -d ' ' -f4)
  srun --dependency=afterok:${job_post} --partition=compute --time=00:15:02 wait.sh

  # run line search
  run_fwat3 $simu_type1
  srun --dependency=afterok:${job_line} --partition=compute --time=00:15:02 wait.sh
  for f in output_fwat3*_${mod}_step*.txt; do mv "$f" "${f/${mod}/${mod}_${simu_type1}}"; done
  run_fwat3 $simu_type2
  srun --dependency=afterok:${job_line} --partition=compute --time=00:15:02 wait.sh
  for f in output_fwat3*_${mod}_step*.txt; do mv "$f" "${f/${mod}/${mod}_${simu_type2}}"; done
  
  # get line search results
  mis1_0=`bash step4_calmisfit.sh 0 $simu_type1 0.0 |cut -d' ' -f1`
  mis2_0=`bash step4_calmisfit.sh 0 $simu_type2 0.0 |cut -d' ' -f1`
  for step in 0.0  `grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`;
  do
    mis1=`bash step4_calmisfit.sh $iter $simu_type1 $step |cut -d' ' -f1`
    mis2=`bash step4_calmisfit.sh $iter $simu_type2 $step |cut -d' ' -f1`
    mis=`echo "$mis1/$mis1_0 + $mis2 / $mis2_0" |bc -l`
    mis=`printf %f $mis`
    echo "$step $mis $mis1 $mis2" >> plots/$mod.mis.avg
  done 

  echo "compute misfit ..."
  fwd=sbash_fwat4_opt_model.sh
  sed -i "/simu_type=/c\simu_type=$simu_type1" $fwd
  sed -i "/MODEL=/c\MODEL=$mod" $fwd
  sed -i "/SETB=/c\SETB=set${setb}" $fwd
  sed -i "/SETE=/c\SETE=set${sete}" $fwd
  job_misfit=$(sbatch $fwd |cut -d ' ' -f4)

  # wait until finish
  srun --dependency=afterok:${job_misfit} --partition=compute --time=00:15:02 wait.sh
  # save log file
  mkdir -p LOGS/$mod
  mv  *.txt LOGS/$mod
done 
