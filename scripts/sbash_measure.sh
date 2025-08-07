#!/bin/bash

# script runs mesher,database generation and solver
# using this example setup 
#
###################################################
source module_env
. parameters.sh

# error flag
set -e 

#==== Comment out the following if running SEM mesh with new models====#
simu_type=noise
NJOBS=8
LOCAL_PC=0

#### STOP HERE #### #
NPROC=`grep ^"NPROC" DATA/Par_file.$simu_type | cut -d'=' -f2`
SOURCE_FILE=src_rec/sources.dat.$simu_type
iter=`fwat-utils getparam iter fwat_params/lbfgs.yaml`
nevts=`awk 'END { print NR }' src_rec/sources.dat.$simu_type`

# mod
MODEL=M`printf %02d $iter`

# working directory
work_dir=`pwd`

# assign job id
TASK_ID=1
if [ "$LOCAL_PC" == "0" ]; then
  TASK_ID=$SLURM_ARRAY_TASK_ID
fi

#logfile
fwd=LOG/output_fwat1_log.$MODEL.$simu_type.job$TASK_ID.txt
FLAG=`fwat-utils getparam flag fwat_params/lbfgs.yaml`
run_opt=3
if [  "$FLAG" == "LS" ]; then 
  run_opt=2
  fwd=LOG/output_fwat3_log.$MODEL.$simu_type.job$TASK_ID.txt
  MODEL=$MODEL.ls
fi
:> $fwd

for i in `seq 1 $NJOBS`; do
  cd $work_dir
  ievt=`echo "($TASK_ID-1) * $NJOBS + $i" |bc`
  ievt_ed=`echo "($TASK_ID-1) * $NJOBS + $NJOBS" |bc`

  # check if job is not included
  if [ "$ievt" -gt  "$nevts" ];then
    break
  fi

  # get evtid
  evtid=`sed -n "$ievt"p $SOURCE_FILE |awk '{print $1}'`
  echo " "
  echo "copying params for $simu_type evtid = $evtid"  

  # prepare files
  fwat-main prepare forward $simu_type $iter $evtid $run_opt

  # run forward simulation
  evtdir=solver/$MODEL/$evtid
  cd $evtdir/
  echo ""
  echo "forward simulation ..."
  date
  $MPIRUN -np $NPROC $SEM_PATH/bin/xspecfem3D
  date

  # merge all seismograms to one big file
  echo "packing seismograms ..."
  fwat-main pack OUTPUT_FILES/seismograms.h5 OUTPUT_FILES/all_seismograms.ascii
  \rm -rf OUTPUT_FILES/all_seismograms.ascii

  # run measure
  echo ""
  echo "measure adjoint source ..."
  cd $work_dir
  date
  $MPIRUN -np $NPROC fwat-main measure $simu_type $iter $evtid $run_opt >> $fwd 
  date

  # adjoint simulation
  fwat-main prepare adjoint $simu_type $iter $evtid $run_opt
  cd $evtdir/
  echo ""
  echo "adjoint simulation ..."
  date
  $MPIRUN -np $NPROC $SEM_PATH/bin/xspecfem3D
  date

  cd $work_dir
  mkdir -p $evtdir/GRADIENT
  \rm -rf $evtdir/GRADIENT/*
  mv $evtdir/DATABASES_MPI/*_kernel.bin $evtdir/GRADIENT
  grad_list=`fwat-model name grad`
  for grad in $grad_list hess_kernel;
  do
    echo "combine $NPROC $grad to hdf5 ..." 
    fwat-main bin2h5 $evtdir/GRADIENT $grad $NPROC 1
  done 
  \rm $evtdir/GRADIENT/*.bin

  # delete useless information
  fwat-utils clean $MODEL $evtid 

  # print flags
  echo " " >> $fwd 
  echo "******************************************************" >> $fwd
  echo "finish event $ievt of $nevts, pair $ievt - $ievt_ed" >> $fwd 
  echo " " >> $fwd
done

