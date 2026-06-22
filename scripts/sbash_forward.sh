###################################################
# error flag
set -e 

source config.env
source utils.sh

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 sbash_forward simu_type"
  exit 1
fi

# input 
simu_type=$1

# params
NPROC=`grep ^"NPROC" DATA/Par_file.$simu_type | cut -d'=' -f2`
SOURCE_FILE=${FWAT_SRC_REC}/sources.dat.$simu_type
iter=`fwat-utils getparam iter ${LBFGS_FILE}`
nevts=`awk 'END { print NR }' ${FWAT_SRC_REC}/sources.dat.$simu_type`
SIMU_TYPES=(`fwat-utils getparam simulation/types| tr -d '[]",'\'`)

# mod
MODEL=M`printf %02d $iter`

# working directory
work_dir=`pwd`

# assign job id
TASK_ID=1
if [ "$PLATFORM"  == "slurm" ];  then 
  TASK_ID=$SLURM_ARRAY_TASK_ID
  for i in "${!SIMU_TYPES[@]}"; do
    if [[ "${SIMU_TYPES[$i]}" == "$simu_type" ]]; then
      NJOBS=${NJOBS_PER_JOBARRAY[$i]}
      break
    fi
  done
elif [ "$PLATFORM"  == "local" ]; then
  TASK_ID=1
  NJOBS=$nevts
else 
  echo "not implemented!"
  exit 1
fi

#logfile
fwd=LOG/output_fwat0_log.$MODEL.$simu_type.job$TASK_ID.txt
:> $fwd

# working dir is IO_TMPDIR when enabled and available
MYDIR=$work_dir
if [ "$USE_IO_TMPDIR" == "1" ];  then
  if [ -d "$IO_TMPDIR" ]; then
    if [ "$PLATFORM"  == "slurm" ];  then
      MYDIR=$IO_TMPDIR/$SLURM_JOB_ID
    elif [ "$PLATFORM"  == "pbs" ];  then
      MYDIR=$IO_TMPDIR/$PBS_JOBID
    fi
    mkdir -p $MYDIR
  else
    USE_IO_TMPDIR=0
  fi
fi
echo "working directory is $MYDIR"
echo "working directory is $MYDIR" >> $fwd

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

  # prepare files
  run_opt=1
  fwat-main prepare forward $simu_type $iter $evtid $run_opt

  # run forward simulation
  evtlist=`cat LOG/.$simu_type-$iter-$evtid-$run_opt`
  \rm LOG/.$simu_type-$iter-$evtid-$run_opt

  # copy fwat_params to MYDIR if using IO_TMPDIR
  if [ "$USE_IO_TMPDIR" == "1" ]; then
    cp -r $work_dir/fwat_params $MYDIR/fwat_params

    # if MOVE_DATABASE is enabled and i == 1, copy model to MYDIR
    if [ "$i" -eq 1 ] && [ "$MOVE_DATABASE" == "1" ]; then
      cp -r $work_dir/$FWAT_OPT_DIR/MODEL_${MODEL} $MYDIR/$MODEL
    fi
  fi

  # run forward simulation
  for evtid_wk in $evtlist;
  do
    evtdir=${FWAT_SOLVER}/$MODEL/$evtid_wk
    if [ "$USE_IO_TMPDIR" == "1" ]; then
      mkdir -p $MYDIR/$evtdir
      cp -r $evtdir/* $MYDIR/$evtdir/

      # remove DATABASES_MPI soft link and link $MYDIR/$MODEL instead
      if [ "$MOVE_DATABASE" == "1" ]; then
        \rm -rf $MYDIR/$evtdir/DATABASES_MPI
        mkdir -p $MYDIR/$evtdir/DATABASES_MPI
        ln -s $MYDIR/${MODEL}/* $MYDIR/$evtdir/DATABASES_MPI/

        # copy axisem data
        if [ -d "$work_dir/DATA/axisem/$evtid_wk" ]; then
          cp -r $work_dir/DATA/axisem/$evtid_wk/* $MYDIR/$evtdir/DATABASES_MPI/
        fi
      fi
    fi

    # go to working directory
    cd $MYDIR/$evtdir/
    echo ""
    echo "forward simulation for $evtid_wk ..."
    date
    $MPIRUN -np $NPROC $SEM_PATH/bin/xspecfem3D
    date

    # merge all seismograms to one big file
    echo "packing seismograms for $evtid_wk ..."
    fwat-main pack OUTPUT_FILES/seismograms.h5 OUTPUT_FILES/all_seismograms.*
    \rm -rf OUTPUT_FILES/all_seismograms.*
    if [ "$USE_IO_TMPDIR" == "1" ]; then
      \cp -r OUTPUT_FILES/seismograms.h5 $work_dir/$evtdir/OUTPUT_FILES/
    fi
    cd $work_dir
  done

  # run measure
  echo ""
  echo "saving forward seismograms for $evtid ..."
  cd $work_dir
  date
  $MPIRUN -np $NPROC_MEASURE fwat-main measure $simu_type $iter $evtid 1 >> $fwd 
  date

  # delete useless information
  for evtid_wk in $evtlist;
  do
    if [ "$USE_IO_TMPDIR" == "1" ]; then
      cd $MYDIR
      fwat-utils clean $MODEL $evtid_wk
    fi
    cd $work_dir
    fwat-utils clean $MODEL $evtid_wk
  done

  # print flags
  echo " " >> $fwd
  echo "******************************************************" >> $fwd
  echo "finish event $ievt of $nevts, pair $ievt - $ievt_ed" >> $fwd
  echo " " >> $fwd
done

# remove
if [ "$USE_IO_TMPDIR" == "1" ]; then
  \rm -rf $MYDIR
fi

