###################################################
# error flag
set -e 

source module_env
. parameters.sh

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 sbash_measure simu_type"
  exit 1
fi


# input 
simu_type=$1

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

  # prepare files
  fwat-main prepare forward $simu_type $iter $evtid $run_opt

  # run forward simulation
  evtlist=`cat LOG/.$simu_type-$iter-$evtid-$run_opt`
  \rm LOG/.$simu_type-$iter-$evtid-$run_opt

  # run forward simulation
  for evtid_wk in $evtlist;
  do 
    evtdir=solver/$MODEL/$evtid_wk
    cd $evtdir/
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
    cd $work_dir
  done 

  # run measure
  echo ""
  echo "measure adjoint source for $evtid ..."
  cd $work_dir
  date
  $MPIRUN -np $NPROC fwat-main measure $simu_type $iter $evtid $run_opt >> $fwd 
  date

  # adjoint simulation
  fwat-main prepare adjoint $simu_type $iter $evtid $run_opt
  for evtid_wk in $evtlist;
  do 
    evtdir=solver/$MODEL/$evtid_wk
    cd $evtdir/
    echo ""
    echo "adjoint simulation for $evtid_wk ..."
    date
    $MPIRUN -np $NPROC $SEM_PATH/bin/xspecfem3D
    date
    echo " "
    cd $work_dir

    # combine kernels
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
    fwat-utils clean $MODEL $evtid_wk 
  done

  # print flags
  echo " " >> $fwd 
  echo "******************************************************" >> $fwd
  echo "finish event $ievt of $nevts, pair $ievt - $ievt_ed" >> $fwd 
  echo " " >> $fwd
done

