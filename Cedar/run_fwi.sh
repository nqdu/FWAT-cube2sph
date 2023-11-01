#!/bin/bash
set -e 

# parameters
simu_type=noise
setb=1
sete=20
iter_start=3
iter_end=6

#simu_type=tele
#setb=1
#sete=22
#iter_start=19
#iter_end=20

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
  echo "iteration $iter"

  # wait util completed
  #if [ "$iter" != "$iter_start" ]; then 
  #  srun --dependency=afterok:${job_misfit} wait.sh
  #fi

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

  # change parfile
  if [ ${mod} == "M00"  ]; then
    LOCAL_PATH="./OUTPUT_FILES/DATABASES_MPI"
    ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH DATA/Par_file
    ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH  DATA/meshfem3D_files/Mesh_Par_file
  else
    LOCAL_PATH="./optimize/MODEL_${mod}"
    ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH DATA/Par_file
    ./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH  DATA/meshfem3D_files/Mesh_Par_file
  fi

  # substitute 
  fwd=sbash_fwat1_fwd_measure_adj.sh
  sed -i "/#SBATCH --array=/c\#SBATCH --array=$setb-$sete" $fwd
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd

  # run forward/adjoint simulation
  echo "computing forward/adjoint simulation ..."
  job_adj=$(sbatch $fwd |cut -d ' ' -f4 )
  #exit

  # run postprocessing
  fwd=sbash_fwat2_postproc_opt.sh 
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/SETB=/c\SETB=set${setb}" $fwd
  sed -i "/SETE=/c\SETE=set${sete}" $fwd 
  echo "postprocessing ..."
  job_post=$(sbatch --dependency=afterok:${job_adj} $fwd | cut -d ' ' -f4)

  # run line search
  if [ $simu_type == "tele" ];then 
    \cp src_rec/sources_ls.dat.tele src_rec/sources_ls.dat -r 
  elif [ $simu_type == "noise"  ]; then
    \cp src_rec/sources_ls.dat.noise src_rec/sources_ls.dat -r 
  else
    \cp src_rec/sources_ls.dat.rf src_rec/sources_ls.dat -r 
  fi
  fwd=sbash_fwat3_linesearch.sh 
  nsteps=`grep NUM_STEP fwat_params/FWAT.PAR |awk -F: '{print $2-1}'`
  sed -i "/#SBATCH --array=/c\#SBATCH --array=0-$nsteps" $fwd
  sed -i "/SIMU_TYPE=/c\SIMU_TYPE=$simu_type" $fwd
  sed -i "/model=/c\model=$mod" $fwd
  echo "line search ..."
  job_line=$(sbatch --dependency=afterok:${job_post} $fwd |cut -d ' ' -f4)
  #exit 
  # compute misfit 
  # compute misfit 
  echo "compute misfit ..."
  fwd=sbash_fwat4_opt_model.sh
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd
  sed -i "/MODEL=/c\MODEL=$mod" $fwd
  sed -i "/SETB=/c\SETB=set${setb}" $fwd
  sed -i "/SETE=/c\SETE=set${sete}" $fwd
  job_misfit=$(sbatch --dependency=afterok:${job_line} $fwd |cut -d ' ' -f4)

  # wait until finish
  srun --dependency=afterok:${job_misfit} --account=rrg-liuqy --time=00:06:02 wait.sh
done 
