#!/bin/bash
set -e 

# parameters
simu_type=rf
setb=23
sete=28
iter_start=0
iter_end=10

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
  fwd=pbs_fwat1_fwd_measure_adj.sh
  sed -i "/#PBS -J/c\#PBS -J $setb-$sete" $fwd
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd

  # run forward/adjoint simulation
  echo "computing forward/adjoint simulation ..."
  job_adj=$(qsub $fwd)

  # run postprocessing
  fwd=pbs_fwat2_postproc_opt.sh 
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/SETB=/c\SETB=set${setb}" $fwd
  sed -i "/SETE=/c\SETE=set${sete}" $fwd 
  echo "postprocessing ..."
  job_post=$(qsub -W depend=afterok:${job_adj} $fwd)

  # run line search
  if [ $simu_type == "tele" ];then 
    \cp src_rec/sources_ls.dat.tele src_rec/sources_ls.dat -r 
  elif [ $simu_type == "noise"  ]; then
    \cp src_rec/sources_ls.dat.noise src_rec/sources_ls.dat -r 
  else
    \cp src_rec/sources_ls.dat.rf src_rec/sources_ls.dat -r 
  fi
  fwd=pbs_fwat3_linesearch.sh 
  nsteps=`grep NUM_STEP fwat_params/FWAT.PAR |awk -F: '{print $2-1}'`
  sed -i "/#PBS -J/c\#PBS -J 0-$nsteps" $fwd
  sed -i "/SIMU_TYPE=/c\SIMU_TYPE=$simu_type" $fwd
  sed -i "/model=/c\model=$mod" $fwd
  echo "line search ..."
  job_line=$(qsub -W depend=afterok:${job_post} $fwd)
  
  # compute misfit 
  echo "compute misfit ..."
  fwd=pbs_fwat4_opt_model.sh
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd
  sed -i "/MODEL=/c\MODEL=$mod" $fwd
  sed -i "/SETB=/c\SETB=set${setb}" $fwd
  sed -i "/SETE=/c\SETE=set${sete}" $fwd
  job_misfit=$(qsub -W depend=afterok:${job_line} $fwd)

  flag=`qstat -x $job_misfit|tail -1 |awk '{print $5}'`
  while [ "$flag" != "F"  ];
  do
     sleep 30
     flag=`qstat -x $job_misfit|tail -1 |awk '{print $5}'`
  done
  mkdir -p LOGS/$mod
  mv  *.o[0-9]* *.txt LOGS/$mod
done 
