###################################################
# error flag
set -e 

source config.env
source utils.sh

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 sbash_measure simu_type"
  exit 1
fi


# input 
simu_type=$1

#### STOP HERE #### #
NPROC=`grep ^"NPROC" DATA/Par_file.$simu_type | cut -d'=' -f2`
SOURCE_FILE=${FWAT_SRC_REC}/sources.dat.$simu_type
iter=`fwat-utils getparam iter ${LBFGS_FILE}`
nevts=`awk 'END { print NR }' ${FWAT_SRC_REC}/sources.dat.$simu_type`
SIMU_TYPES=(`fwat-utils getparam simulation/types| tr -d '[]",'\'`)

# mod
MODEL=M`printf %02d $iter`

# current directory
curr_dir=`pwd`

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
fwd=LOG/output_fwat1_log.$MODEL.$simu_type.job$TASK_ID.txt
FLAG=`fwat-utils getparam flag ${LBFGS_FILE}`
run_opt=3
if [  "$FLAG" == "LS" ]; then 
  run_opt=2
  fwd=LOG/output_fwat3_log.$MODEL.$simu_type.job$TASK_ID.txt
  MODEL=$MODEL.ls
fi
:> $fwd
echo " " >> $fwd

# working dir is IO_TMPDIR when enabled and available
# command to run one-node command on compute nodes
# e.g. pbsdsh -u for PBS
# e.g. srun --ntasks-per-node=1 for slurm
# e.g. mpirun --map-by ppr:1:node for slurm/pbs with mpirun, works on Scinet
# e.g empty for local
MYDIR=$curr_dir
PRUN=""
if [ "$USE_IO_TMPDIR" == "1" ];  then
  if [ -d "$IO_TMPDIR" ]; then
    if [ "$PLATFORM"  == "slurm" ];  then 
      MYDIR=$IO_TMPDIR/$SLURM_JOB_ID
    elif [ "$PLATFORM"  == "pbs" ];  then 
      MYDIR=$IO_TMPDIR/$PBS_JOBID
    fi
    PRUN="mpirun --map-by ppr:1:node" 
    $PRUN mkdir -p $MYDIR
  else
    USE_IO_TMPDIR=0
  fi
fi
echo "working directory is $MYDIR" 
echo "working directory is $MYDIR" >> $fwd

for i in `seq 1 $NJOBS`; do
  cd $curr_dir
  ievt=`echo "($TASK_ID-1) * $NJOBS + $i" |bc`
  ievt_ed=`echo "($TASK_ID-1) * $NJOBS + $NJOBS" |bc`

  # check if job is not included
  if [ "$ievt" -gt  "$nevts" ];then
    break
  fi

  # get evtid
  evtid=`sed -n "$ievt"p $SOURCE_FILE |awk '{print $1}'`

  # prepare files
  fwat-main prepare forward $simu_type $iter $evtid $run_opt

  # run forward simulation
  evtlist=`cat LOG/.$simu_type-$iter-$evtid-$run_opt`
  \rm LOG/.$simu_type-$iter-$evtid-$run_opt

  # copy fwat_params to MYDIR if using IO_TMPDIR
  if [ "$USE_IO_TMPDIR" == "1" ]; then
    $PRUN \cp -r $curr_dir/fwat_params $MYDIR/fwat_params

    # if MOVE_DATABASE is enabled and i == 1, copy model to MYDIR
    if [ "$i" -eq 1 ] && [ "$MOVE_DATABASE" == "1" ]; then
      $PRUN \cp -r $curr_dir/$FWAT_OPT_DIR/MODEL_${MODEL} $MYDIR/$MODEL
    fi
  fi

  # run forward simulation
  for evtid_wk in $evtlist;
  do 
    evtdir=${FWAT_SOLVER}/$MODEL/$evtid_wk
    if [ "$USE_IO_TMPDIR" == "1" ]; then
      $PRUN mkdir -p $MYDIR/$evtdir
      $PRUN \cp -r $evtdir/* $MYDIR/$evtdir/

      # remove DATABASES_MPI soft link and link $MYDIR/$MODEL instead
      if [ "$MOVE_DATABASE" == "1" ]; then
        \rm -rf $MYDIR/$evtdir/DATABASES_MPI
        mkdir -p $MYDIR/$evtdir/DATABASES_MPI
        ln -s $MYDIR/${MODEL}/* $MYDIR/$evtdir/DATABASES_MPI/

        # copy axisem data
        if [ -d "$curr_dir/DATA/axisem/$evtid_wk" ]; then
          $PRUN \cp -r $curr_dir/DATA/axisem/$evtid_wk/* $MYDIR/$evtdir/DATABASES_MPI/
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
    echo " "
    echo "packing seismograms for $evtid_wk ..."
    fwat-main pack OUTPUT_FILES/seismograms.h5 OUTPUT_FILES/all_seismograms.*
    \rm -rf OUTPUT_FILES/all_seismograms.*
    if [ "$USE_IO_TMPDIR" == "1" ]; then
      # copy back to solver dir, as seismograms.h5 is on primary compute node.
      \cp -r OUTPUT_FILES/seismograms.h5 $curr_dir/$evtdir/OUTPUT_FILES/
    fi
    \cp $MYDIR/$evtdir/OUTPUT_FILES/output_solver.txt $curr_dir/$evtdir/OUTPUT_FILES/output_solver.fwd.txt

    # go back to current dir 
    cd $curr_dir
  done 

  # run measure
  echo ""
  echo "measure adjoint source for $evtid ..."
  cd $curr_dir
  date
  $MPIRUN -np $NPROC_MEASURE fwat-main measure $simu_type $iter $evtid $run_opt >> $fwd 
  date

  # adjoint simulation
  fwat-main prepare adjoint $simu_type $iter $evtid $run_opt
  for evtid_wk in $evtlist;
  do 
    evtdir=${FWAT_SOLVER}/$MODEL/$evtid_wk
    if [ "$USE_IO_TMPDIR" == "1" ]; then
      \rm -rf $MYDIR/$evtdir/SEM
      $PRUN cp -r $curr_dir/$evtdir/SEM $MYDIR/$evtdir/
      \rm -rf $curr_dir/$evtdir/SEM 
      $PRUN \cp -r $curr_dir/$evtdir/DATA/* $MYDIR/$evtdir/DATA/
      cd $MYDIR/$evtdir
    fi

    # go to working directory
    cd $MYDIR/$evtdir
    echo ""
    echo "adjoint simulation for $evtid_wk ..."
    date
    $MPIRUN -np $NPROC $SEM_PATH/bin/xspecfem3D
    date
    echo " "
    
    # combine kernels, note we are on NLS
    cd $MYDIR
    $PRUN mkdir -p $evtdir/GRADIENT
    $PRUN bash -c "find $evtdir/GRADIENT/ -maxdepth 1 -name 'proc*_kernel.bin' -print0 | xargs -0 rm -f"
    find $evtdir/DATABASES_MPI/ -maxdepth 1 -name 'proc*_kernel.bin' -print0 | xargs -0 mv -t $evtdir/GRADIENT/
    mpirun -np $NPROC fwat-model combine_kl $evtdir/DATABASES_MPI/ $evtdir/GRADIENT
    find $evtdir/GRADIENT/ -maxdepth 1 -name 'proc*_kernel.bin' -print0 | xargs -0 rm -f
    rm -rf $curr_dir/$evtdir/GRADIENT
    \cp -r $MYDIR/$evtdir/GRADIENT $curr_dir/$evtdir/
    $PRUN \rm -rf $MYDIR/$evtdir/GRADIENT
    echo ""

    # delete useless information
    if [ "$USE_IO_TMPDIR" == "1" ]; then
      $PRUN fwat-utils clean $MODEL $evtid_wk 
    fi

    # copy back to solver dir
    cd $curr_dir
    fwat-utils clean $MODEL $evtid_wk 
    \cp $MYDIR/$evtdir/OUTPUT_FILES/output_solver.txt $curr_dir/$evtdir/OUTPUT_FILES/output_solver.adj.txt
  done

  # print flags
  echo " " >> $fwd 
  echo "******************************************************" >> $fwd
  echo "finish event $ievt of $nevts, pair $ievt - $ievt_ed" >> $fwd 
  echo " " >> $fwd
done

# remove
if [ "$USE_IO_TMPDIR" == "1" ]; then
  $PRUN \rm -rf $MYDIR
fi

