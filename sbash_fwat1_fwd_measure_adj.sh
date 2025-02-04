#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=160
#SBATCH --array=21-28%5
#SBATCH --time=00:35:59
#SBATCH --job-name FWD_ADJ
#SBATCH --output=FWD_ADJ-%j_set%a.txt
#SBATCH --partition=compute
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
source module_env
export OMP_NUM_THREADS=1

# error flag
set -e 

# include file
. parameters.sh

#==== Comment out the following if running SEM mesh with new models====#
MODEL=M00
simu_type=tele
NJOBS=4
START_SET=1
NPROC=$SLURM_NTASKS

# parfile changer script
change_par=$FWATLIB/change_par_file.sh

SOURCE_FILE=src_rec/sources.dat.$simu_type
iter=`echo $MODEL |cut -d"M" -f2 |awk '{printf "%d", $1}'`
nevts=`cat $SOURCE_FILE |wc -l`
work_dir=`pwd`
mod=$MODEL

#logfile
fwd=output_fwat1_log.$MODEL.$simu_type.job$SLURM_ARRAY_TASK_ID.txt
:> $fwd

for i in `seq 1 $NJOBS`; do
  cd $work_dir
  ievt=`echo "($SLURM_ARRAY_TASK_ID-1) * $NJOBS + $i" |bc`
  ievt_ed=`echo "($SLURM_ARRAY_TASK_ID-1) * $NJOBS + $NJOBS" |bc`
  id=`echo "$START_SET + $ievt -1" |bc`

  # check if job is not included
  if [ "$ievt" -gt  "$nevts" ];then
    break
  fi

  # get evtid
  evtid=`sed -n "$ievt"p $SOURCE_FILE |awk '{print $1}'`
  echo "copying params for $simu_type evtid = $evtid"  
  evtdir=$work_dir/solver/$mod/$evtid

  # mkdir
  mkdir -p $evtdir/
  cd $evtdir/
  mkdir -p DATA  OUTPUT_FILES/  DATABASES_MPI DATA/meshfem3D_files
  cd $work_dir
  
  # copy common parameters
  \cp DATA/Par_file.$simu_type $evtdir/DATA/Par_file
  \cp fwat_params/FWAT.PAR.$simu_type $evtdir/DATA/FWAT.PAR 
  \cp fwat_params/MEASUREMENT.PAR.$simu_type $evtdir/DATA/MEASUREMENT.PAR
  \cp OUTPUT_FILES/*.h $evtdir/OUTPUT_FILES
  \rm -rf $evtdir/DATA/meshfem3D_files/*
  ln -s $work_dir/DATA/meshfem3D_files/* $evtdir/DATA/meshfem3D_files/
  \cp DATA/adepml_stage $evtdir/DATA/
  \cp DATA/wavefield* $evtdir/DATA/

  # model link
  LOCAL_PATH="./DATABASES_MPI"
  $change_par LOCAL_PATH $LOCAL_PATH $evtdir/DATA/Par_file
  $change_par LOCAL_PATH $LOCAL_PATH  $evtdir/DATA/meshfem3D_files/Mesh_Par_file
  \rm -rf $evtdir/$LOCAL_PATH/*
  if [ ${mod} == "M00"  ]; then
    ln -s $work_dir/$LOCAL_PATH/* $evtdir/$LOCAL_PATH/
  else
    ln -s $work_dir/./optimize/MODEL_${mod}/* $evtdir/$LOCAL_PATH/
  fi

  # copy simu_type depedent files
  cd src_rec
  \cp STATIONS_$evtid $evtdir/DATA/STATIONS
  \cp STATIONS_$evtid $evtdir/DATA/STATIONS_ADJOINT
  if [[ $simu_type == "tele" ]]; then 
    :> $evtdir/DATA/FORCESOLUTION
    
    # change parameters
    cd $work_dir
    $change_par COUPLE_WITH_INJECTION_TECHNIQUE .true. $evtdir/DATA/Par_file
    $change_par INJECTION_TECHNIQUE_TYPE 4 $evtdir/DATA/Par_file

    # link axisem field
    \rm -rf $evtdir/$LOCAL_PATH/*wavefield_discontinuity.bin
    ln -s $work_dir/DATA/axisem/$evtid/* $evtdir/$LOCAL_PATH/
  elif [[ $simu_type == "noise" ]]; then
    \cp FORCESOLUTION_$evtid $evtdir/DATA/FORCESOLUTION

    cd ..
    $change_par COUPLE_WITH_INJECTION_TECHNIQUE .false. $evtdir/DATA/Par_file
  else 
    echo "noise"
  fi
  cd $work_dir

  # forward simulation
  echo ""
  echo "forward simulation ..."
  $change_par SUBSAMPLE_FORWARD_WAVEFIELD .true. $evtdir/DATA/Par_file
  $change_par SIMULATION_TYPE 1 $evtdir/DATA/Par_file
  $change_par APPROXIMATE_HESS_KL .false. $evtdir/DATA/Par_file
  cd $evtdir/
  mpirun -np $NPROC $fksem/bin/xspecfem3D

  # run measure
  echo ""
  cd $work_dir
  bash $MEASURE_LIB/measure.$simu_type.sh $iter $evtid 3 >> $fwd 

  # adjoint simulation
  $change_par COUPLE_WITH_INJECTION_TECHNIQUE .false. $evtdir/DATA/Par_file
  $change_par SUBSAMPLE_FORWARD_WAVEFIELD .true. $evtdir/DATA/Par_file
  $change_par SIMULATION_TYPE 3 $evtdir/DATA/Par_file
  $change_par APPROXIMATE_HESS_KL .true. $evtdir/DATA/Par_file
  cd $evtdir/
  echo ""
  echo "adjoint simulation ..."
  mpirun -np $NPROC $fksem/bin/xspecfem3D

  # copy kernels to GRADIENT
  cd $work_dir
  mkdir -p $evtdir/GRADIENT
  \rm -rf $evtdir/GRADIENT/*
  mv $evtdir/$LOCAL_PATH/*_kernel.bin $evtdir/GRADIENT

  # delete useless information
  bash $FWATLIB/clean.sh $mod $evtid 

  # print flags
  echo " " >> $fwd 
  echo "******************************************************" >> $fwd
  echo "finish event $ievt of $nevts, pair $ievt - $ievt_ed" >> $fwd 
  echo " " >> $fwd
done

