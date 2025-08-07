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
  echo " "
  echo "copying params for $simu_type evtid = $evtid"  

  # prepare files
  fwat-main prepare forward $simu_type $iter $evtid 1

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
  echo "saving forward seismograms ..."
  cd $work_dir
  date
  $MPIRUN -np $NPROC fwat-main measure $simu_type $iter $evtid 1 >> $fwd 
  date

  # delete useless information
  fwat-utils clean $MODEL $evtid 

  # print flags
  echo " " >> $fwd 
  echo "******************************************************" >> $fwd
  echo "finish event $ievt of $nevts, pair $ievt - $ievt_ed" >> $fwd 
  echo " " >> $fwd
done

