#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=1-4%5
#SBATCH --time=00:35:59
#SBATCH --job-name=FWD
#SBATCH --output=FWD-%j_set%a.txt
#SBATCH --account=rrg-liuqy
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca

# script runs mesher,database generation and solver
# using this example setup 
#
###################################################
source module_env
. parameters.sh

# error flag
set -e 

#==== Comment out the following if running SEM mesh with new models====#
MODEL=M00
simu_type=noise
NJOBS=8
START_SET=1
LOCAL_PC=0

#### STOP HERE #### #

NPROC=`grep ^"NPROC" DATA/Par_file.$simu_type | cut -d'=' -f2`


# parfile changer script
change_par=$FWATLIB/change_par_file.sh

SOURCE_FILE=src_rec/sources.dat.$simu_type
iter=`echo $MODEL |cut -d"M" -f2 |awk '{printf "%d", $1}'`
nevts=`awk 'END { print NR }' src_rec/sources.dat.$simu_type`
work_dir=`pwd`
mod=$MODEL

# assign job id
TASK_ID=1
if [ "$LOCAL_PC" == "0" ]; then
  TASK_ID=$SLURM_ARRAY_TASK_ID
fi

#logfile
fwd=output_fwat0_log.$MODEL.$simu_type.job$TASK_ID.txt
:> $fwd

for i in `seq 1 $NJOBS`; do
  cd $work_dir
  ievt=`echo "($TASK_ID-1) * $NJOBS + $i" |bc`
  ievt_ed=`echo "($TASK_ID-1) * $NJOBS + $NJOBS" |bc`
  id=`echo "$START_SET + $ievt -1" |bc`

  # check if job is not included
  if [ "$ievt" -gt  "$nevts" ];then
    break
  fi

  # get evtid
  evtid=`sed -n "$ievt"p $SOURCE_FILE |awk '{print $1}'`
  echo " "
  echo "copying params for $simu_type evtid = $evtid"  
  evtdir=$work_dir/solver/$mod/$evtid

  # mkdir
  mkdir -p $evtdir/
  cd $evtdir/
  mkdir -p DATA  OUTPUT_FILES/  DATABASES_MPI DATA/meshfem3D_files
  cd $work_dir
  
  # copy common parameters
  \cp DATA/Par_file.$simu_type $evtdir/DATA/Par_file
  \cp fwat_params/FWAT.PAR.yaml $evtdir/DATA/FWAT.PAR.yaml
  \cp OUTPUT_FILES/*.h $evtdir/OUTPUT_FILES
  \rm -rf $evtdir/DATA/meshfem3D_files/*
  ln -s $work_dir/DATA/meshfem3D_files/* $evtdir/DATA/meshfem3D_files/

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
  if [[ $simu_type == "tele" ]] || \
     [[ $simu_type == "sks" ]] ; then 
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
    echo "not implemented!"
    exit 1
  fi
  cd $work_dir

  # forward simulation
  echo ""
  echo "forward simulation ..."
  $change_par SUBSAMPLE_FORWARD_WAVEFIELD .false. $evtdir/DATA/Par_file
  $change_par SIMULATION_TYPE 1 $evtdir/DATA/Par_file
  $change_par APPROXIMATE_HESS_KL .false. $evtdir/DATA/Par_file
  cd $evtdir/
  date
  $MPIRUN -np $NPROC $fksem/bin/xspecfem3D
  date

  # merge all seismograms to one big file
  echo "packing seismograms ..."
  python $MEASURE_LIB/pack_seismogram.py OUTPUT_FILES/seismograms.h5 OUTPUT_FILES/*.semd
  \rm -rf OUTPUT_FILES/*.semd

  # run measure
  echo ""
  echo "saving forward seismograms ..."
  cd $work_dir
  date
  $MPIRUN -np $NPROC python $MEASURE_LIB/run_preprocess.py $simu_type $iter $evtid 1 >> $fwd 
  date

  # delete useless information
  bash $FWATLIB/clean.sh $mod $evtid 

  # print flags
  echo " " >> $fwd 
  echo "******************************************************" >> $fwd
  echo "finish event $ievt of $nevts, pair $ievt - $ievt_ed" >> $fwd 
  echo " " >> $fwd
done

