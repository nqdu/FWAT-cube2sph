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

# working dir is IO_TMPDIR when enabled and available
MYDIR=$curr_dir
if [ "$USE_IO_TMPDIR" == "1" ];  then
  if [ -d "$IO_TMPDIR" ]; then
    echo "working directory is $IO_TMPDIR"
    if [ "$PLATFORM"  == "slurm" ];  then 
      MYDIR=$IO_TMPDIR/$SLURM_JOB_ID
    elif [ "$PLATFORM"  == "pbs" ];  then 
      MYDIR=$IO_TMPDIR/$PBS_JOBID
    fi
    mkdir -p $MYDIR
  else
    echo "IO_TMPDIR is not available, disable USE_IO_TMPDIR"
    USE_IO_TMPDIR=0
    echo "working directory is current dir"
  fi
fi

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
    cp -r $curr_dir/fwat_params $MYDIR/fwat_params
  fi

  # run forward simulation
  for evtid_wk in $evtlist;
  do 
    evtdir=${FWAT_SOLVER}/$MODEL/$evtid_wk
    if [ "$USE_IO_TMPDIR" == "1" ]; then
      mkdir -p $MYDIR/$evtdir
      cp -r $evtdir/* $MYDIR/$evtdir/
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
      mv $curr_dir/$evtdir/SEM $MYDIR/$evtdir/
      \cp -r $curr_dir/$evtdir/DATA/* $MYDIR/$evtdir/DATA/
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

    # combine kernels
    cd $MYDIR
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
    echo ""

    # delete useless information
    if [ "$USE_IO_TMPDIR" == "1" ]; then
      fwat-utils clean $MODEL $evtid_wk 
    fi

    # copy back to solver dir
    cd $curr_dir
    fwat-utils clean $MODEL $evtid_wk 
    if [ "$USE_IO_TMPDIR" == "1" ]; then
      \rm -rf $evtdir/GRADIENT/
      mv $MYDIR/$evtdir/GRADIENT $evtdir/
    fi
    \cp $MYDIR/$evtdir/OUTPUT_FILES/output_solver.txt $curr_dir/$evtdir/OUTPUT_FILES/output_solver.adj.txt
  done

  # print flags
  echo " " >> $fwd 
  echo "******************************************************" >> $fwd
  echo "finish event $ievt of $nevts, pair $ievt - $ievt_ed" >> $fwd 
  echo " " >> $fwd
done

