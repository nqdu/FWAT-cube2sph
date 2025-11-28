###################################################
# error flag
set -e 

source module_env
. parameters.sh

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

  # run forward simulation
  nsta_used=0
  for evtid_wk in $evtlist;
  do 
    evtdir=${FWAT_SOLVER}/$MODEL/$evtid_wk
    cd $evtdir/
    echo ""
    echo "forward simulation for $evtid_wk ..."
    date
    $MPIRUN -np $NPROC $SEM_PATH/bin/xspecfem3D
    date

    # merge all seismograms to one big file
    echo "packing seismograms for $evtid_wk ..."
    fwat-main pack OUTPUT_FILES/seismograms.h5 OUTPUT_FILES/all_seismograms.*
    \rm -rf OUTPUT_FILES/all_seismograms.*
    cd $work_dir

    # count stations used
    nsta_used=`awk 'END{print NR}' ${FWAT_SRC_REC}/DATA/STATIONS`
  done 

  # run measure
  echo ""
  echo "saving forward seismograms for $evtid ..."
  cd $work_dir
  date
  local nproc_run=$NPROC 
  if [ $nsta_used -lt $NPROC ]; then
    nproc_run=$nsta_used
  fi
  $MPIRUN -np $nproc_run fwat-main measure $simu_type $iter $evtid 1 >> $fwd 
  date

  # delete useless information
  for evtid_wk in $evtlist;
  do
    fwat-utils clean $MODEL $evtid_wk 
  done  

  # print flags
  echo " " >> $fwd 
  echo "******************************************************" >> $fwd
  echo "finish event $ievt of $nevts, pair $ievt - $ievt_ed" >> $fwd 
  echo " " >> $fwd
done

