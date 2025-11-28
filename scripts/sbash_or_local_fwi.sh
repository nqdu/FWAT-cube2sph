#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=00:30:00
#SBATCH --job-name=ADJ
#SBATCH --output=LOG/ADJ-%j.txt
#SBATCH --account=rrg-liuqy
#SBATCH --partition=compute
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca

###################################################
# error flag
set -e 

run_one_simu_() {
  local simu_type=$1
  local evtid=$2
  local iter=$3
  local run_opt=$4 
  local fwd=$5

  # prepare files
  fwat-main prepare forward $simu_type $iter $evtid $run_opt

  # run forward simulation
  local evtlist=`cat LOG/.$simu_type-$iter-$evtid-$run_opt`
  \rm LOG/.$simu_type-$iter-$evtid-$run_opt

  # run forward simulation and count stations used
  local nsta_used=0
  for evtid_wk in $evtlist;
  do 
    evtdir=${FWAT_SOLVER}/$MODEL/$evtid_wk
    cd $evtdir/
    echo ""
    echo "forward simulation for $evtid_wk `date` ..."
    date
    $MPIRUN -np $NPROC $SEM_PATH/bin/xspecfem3D
    echo "finished $evtid_wk at `date`"

    # merge all seismograms to one big file
    echo " "
    echo "packing seismograms for $evtid_wk ..."
    fwat-main pack OUTPUT_FILES/seismograms.h5 OUTPUT_FILES/all_seismograms.*
    \rm -rf OUTPUT_FILES/all_seismograms.*
    cd $work_dir

    # count stations used
    nsta_used=`awk 'END{print NR}' ${FWAT_SRC_REC}/DATA/STATIONS`
  done 

  # run measure
  echo ""
  echo "measure adjoint source for $evtid at `date` ..."
  cd $work_dir
  local nproc_run=$NPROC 
  if [ $nsta_used -lt $NPROC ]; then
    nproc_run=$nsta_used
  fi
  $MPIRUN -np $nproc_run fwat-main measure $simu_type $iter $evtid $run_opt >> $fwd 
  echo "finished measure for $evtid at `date`"

  # adjoint simulation
  fwat-main prepare adjoint $simu_type $iter $evtid $run_opt
  for evtid_wk in $evtlist;
  do 
    evtdir=${FWAT_SOLVER}/$MODEL/$evtid_wk
    cd $evtdir/
    echo ""
    echo "adjoint simulation for $evtid_wk at `date` ..."
    $MPIRUN -np $NPROC $SEM_PATH/bin/xspecfem3D
    echo "finished adjoint for $evtid_wk at `date`"
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

}

run_measure()
{
  local iter=$1
  local NPROCS_TOTAL=$2

  # check copy all events and simutype into arrays
  local evtid_list=()
  local simu_type_list=()
  for simu_type in "${SIMU_TYPES[@]}"; do
    local SOURCE_FILE=${FWAT_SRC_REC}/sources.dat.$simu_type
    evtid_list+=(`awk '{print $1}' $SOURCE_FILE`)
    simu_type_list+=(`awk -v a=$simu_type '{print $1,a}' $SOURCE_FILE | awk '{print $2}'`)
  done

  NPROC=`grep ^"NPROC" DATA/Par_file.$simu_type | cut -d'=' -f2` # NPROCS per simu
  nsimus_parallel=`echo "$NPROCS_TOTAL / $NPROC" | bc`
  if [ $nsimus_parallel -lt 1 ]; then
    # exit if not enough procs
    echo "Error: not enough total procs $NPROCS_TOTAL for one simulation requiring $NPROC procs"
    exit 1
  fi

  # mod
  local FLAG=`fwat-utils getparam flag ${LBFGS_FILE}`
  local MODEL=M`printf %02d $iter`
  local run_opt=3
  if [  "$FLAG" == "LS" ]; then 
    local run_opt=2
    local fwd=LOG/output_fwat3_log.$MODEL.$simu_type.$i.txt
    MODEL=$MODEL.ls
  fi

  # log files
  local fwd_tag=LOG/output_fwat1_log.$MODEL.$simu_type
  for ((i=0; i < $nsimus_parallel; i++)); do
    local simu_type=${SIMU_TYPES[0]}
    
    if [  "$FLAG" == "LS" ]; then 
      fwd_tag=LOG/output_fwat3_log.$MODEL.$simu_type
    fi
    :> $fwd_tag.$i.txt
  done

  # loop over all events 
  local nevts=${#evtid_list[@]}
  for (( i=0; i<$nevts; i+=nsimus_parallel )); do
    # manually launch simu in parallel
    for (( j=0; j<$nsimus_parallel; j++ )); do
      local idx=$(( i + j ))
      if [ "$idx" -ge "$nevts" ]; then
        break
      fi

      # get id
      local evtid=${evtid_list[$idx]}
      local simu_type=${simu_type_list[$idx]}
      #echo "launching event $evtid $simu_type  iter $iter  run_opt $run_opt  in background ..."
      run_one_simu_ $simu_type $evtid $iter $run_opt $fwd_tag.$j.txt &
    done 

    # wait for all to finish
    wait
    
  done 
}

source module_env
. parameters.sh

# submit and get job id
if [ "$PLATFORM"  == "local"  ]; then 
  if [ "$#" -ne 2 ]; then
    echo "Usage: $0 sbash_measure NPROCS_TOTAL max_iter"
    exit 1
  fi
  NPROCS_TOTAL=$1
  max_iter=$2
else 
  if [ "$#" -ne 1 ]; then
    echo "Usage: $0 sbash_measure max_iter"
    exit 1
  fi
  NPROCS_TOTAL=$(( SLURM_NNODES * SLURM_NTASKS_PER_NODE ))
  max_iter=$1
fi

# working directory
work_dir=`pwd`

# create working directories
mkdir -p misfits optimize solver LOG

# now loop over iterations 
for ii in `seq 1 $max_iter`;do 

  # current model
  iter=`fwat-utils getparam iter ${LBFGS_FILE}`
  flag=`fwat-utils getparam flag ${LBFGS_FILE}`
  mod=M`printf %02d $iter`
  echo ""
  echo "iteration $iter $mod $flag"
  echo "============================"
  echo " "
  

  # check flag type and run 
  if [ $flag == "INIT" ]; then 
    run_measure $iter $NPROCS_TOTAL
    #exit 
    
    # sum kernels, get search direction, generate trial model 
    bash sbash_postproc_kl.sh > LOG/POST.$iter.txt

  elif [ $flag == "GRAD"  ];then 
    # get search direction, generate trial model 
    bash sbash_postproc_kl.sh > LOG/POST.$iter.txt

  else  # line search
    run_measure $iter $NPROCS_TOTAL

    bash sbash_wolfe.sh > LOG/WOLFE.$iter.txt
  fi

done


