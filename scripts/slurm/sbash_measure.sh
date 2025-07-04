#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --array=21-28%5
#SBATCH --time=00:35:59
#SBATCH --job-name=FWD_ADJ
#SBATCH --output=FWD_ADJ-%j_set%a.txt
#SBATCH --account=rrg-liuqy
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca

# error flag
set -e 

#
###################################################
source module_env 
. parameters.sh

#==== Comment out the following if running SEM mesh with new models====#
MODEL=M00
simu_type=tele
NJOBS=4
START_SET=1
NPROC=`grep ^"NPROC" DATA/Par_file.$simu_type | cut -d'=' -f2`
LOCAL_PC=1

# parfile changer script
change_par=$FWATLIB/change_par_file.sh

SOURCE_FILE=src_rec/sources.dat.$simu_type
iter=`echo $MODEL |cut -d"M" -f2 |awk '{printf "%d", $1}'`
nevts=`awk 'END { print NR }' src_rec/sources.dat.$simu_type`
mod=$MODEL
work_dir=`pwd`

# assign job id
TASK_ID=1
if [ "$LOCAL_PC" == "0" ]; then
  TASK_ID=$SLURM_ARRAY_TASK_ID
fi

#logfile
fwd=LOG/output_fwat1_log.$MODEL.$simu_type.job$TASK_ID.txt

# check for LS or INIT
FLAG=`python $FWATLIB/get_param.py flag $FWATPARAM/lbfgs.yaml`
run_opt=3
if [  "$FLAG" == "LS" ]; then 
  run_opt=2
  fwd=LOG/output_fwat3_log.$MODEL.$simu_type.job$TASK_ID.txt
  mod=$MODEL.ls
fi
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
  echo "job $i of $NJOBS : copying params for $simu_type evtid = $evtid"  
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
  ln -s $work_dir/./optimize/MODEL_${mod}/* $evtdir/$LOCAL_PATH/

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
  $change_par SUBSAMPLE_FORWARD_WAVEFIELD .true. $evtdir/DATA/Par_file
  $change_par SIMULATION_TYPE 1 $evtdir/DATA/Par_file
  $change_par APPROXIMATE_HESS_KL .false. $evtdir/DATA/Par_file
  $change_par WRITE_SEISMOGRAMS_BY_MASTER .true. $evtdir/DATA/Par_file
  $change_par SAVE_ALL_SEISMOS_IN_ONE_FILE .true. $evtdir/DATA/Par_file
  NSTEP=`grep '^NSTEP ' $evtdir/DATA/Par_file |awk -F'=' '{print $2}'`
  $change_par NTSTEP_BETWEEN_OUTPUT_SEISMOS $NSTEP $evtdir/DATA/Par_file
  cd $evtdir/
  date
  $MPIRUN -np $NPROC $fksem/bin/xspecfem3D
  date

  # merge all seismograms to one big file
  echo "packing seismograms ..."
  python $MEASURE_LIB/pack_seismogram.py OUTPUT_FILES/seismograms.h5 OUTPUT_FILES/all_seismograms.ascii
  \rm -rf OUTPUT_FILES/all_seismograms.ascii

  # run measure
  echo ""
  echo "measure adjoint source ..."
  cd $work_dir
  date
  $MPIRUN -np $NPROC python $MEASURE_LIB/run_preprocess.py $simu_type $iter $evtid $run_opt >> $fwd
  date

  # adjoint simulation
  $change_par COUPLE_WITH_INJECTION_TECHNIQUE .false. $evtdir/DATA/Par_file
  $change_par SUBSAMPLE_FORWARD_WAVEFIELD .true. $evtdir/DATA/Par_file
  $change_par SIMULATION_TYPE 3 $evtdir/DATA/Par_file
  $change_par APPROXIMATE_HESS_KL .true. $evtdir/DATA/Par_file
  cd $evtdir/
  echo ""
  echo "adjoint simulation ..."
  date
  $MPIRUN -np $NPROC $fksem/bin/xspecfem3D
  date

  # copy kernels to GRADIENT
  # and combine them to hdf5
  cd $work_dir
  mkdir -p $evtdir/GRADIENT
  \rm -rf $evtdir/GRADIENT/*
  mv $evtdir/$LOCAL_PATH/*_kernel.bin $evtdir/GRADIENT
  grad_list=`GET_GRAD_NAME`
  for grad in $grad_list hess_kernel;
  do
    echo "combine $NPROC $grad to hdf5 ..." 
    python $OPT_LIB/bin2h5.py $evtdir/GRADIENT $grad $NPROC 1
  done 
  \rm $evtdir/GRADIENT/*.bin

  # delete useless information
  bash $FWATLIB/clean.sh $mod $evtid 

  # print flags
  echo " " >> $fwd 
  echo "******************************************************" >> $fwd
  echo "finish event $ievt of $nevts, pair $ievt - $ievt_ed" >> $fwd 
  echo " " >> $fwd
done

