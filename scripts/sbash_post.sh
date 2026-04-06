#!/bin/bash 
set -e 
run_wolfe () { # run wolfe line search
  # compute misfit
  if [ "$nsimtypes" == "1" ]; then 
    # compute misfits
    info=`fwat-main misfit $MODEL ${SIMU_TYPES[0]}`
    chi=`echo $info |awk '{print $1}'`
    info=`fwat-main misfit $MODEL.ls ${SIMU_TYPES[0]}`
    chi1=`echo $info |awk '{print $1}'`
  else
    #init misifits
    chi=0.
    chi1=0.

    # for all simulation types, compute weighted misfits
    for((i=0;i<$nsimtypes;i++)); 
    do 
      info=`fwat-main misfit $MODEL ${SIMU_TYPES[$i]}`
      l1=`echo $info |awk '{print $1}'`
      info=`fwat-main misfit $MODEL.ls ${SIMU_TYPES[$i]}`
      l2=`echo $info |awk '{print $1}'`

      # weighted sum for iteration n
      # \sum_i L^(n)_i / L^(0)_i * L^(0)_0 * user_weight
      chi=`echo $chi $l1 ${SIMU_WEIGHTS[$i]} |awk '{print $1+$2*$3}'`
      chi1=`echo $chi1 $l2 ${SIMU_WEIGHTS[$i]}|awk '{print $1+$2*$3}'`
    done 
  fi 

  echo "misfit current/next = $chi $chi1"
  echo " "

  logfile=LOG/output_fwat4_log_$MODEL.txt
  echo "******************************************************" > $logfile
  
  # get smooth parameters
  GPU_MODE=`grep ^"GPU_MODE" DATA/Par_file | cut -d'=' -f2`
  LOCAL_PATH=${FWAT_OPT_DIR}/MODEL_${MODEL}
  change_par LOCAL_PATH $LOCAL_PATH ./DATA/Par_file
  change_par LOCAL_PATH $LOCAL_PATH ./DATA/meshfem3D_files/Mesh_Par_file
  info=`fwat-utils getparam optimize/SMOOTHING  | sed 's/\[\|]//g' | sed 's/,/ /g'`
  sigma_h=`echo $info | awk  '{print $1}'`
  sigma_v=`echo $info | awk  '{print $2}'`

  # sum kernels for line search, save to optimize/sum_kernels_$MODEL.ls
  echo "sum kernels for new model ..."
  $MPIRUN $hostfile -np $NPROC fwat-main sum_kernel $MODEL.ls
  echo " " >> $logfile

  # smooth gradient 
  kl_list=`fwat-model name direc`
  for param1 in $kl_list; 
  do 
    param=${param1:1}_kernel
    mv ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}.ls/*_${param}.bin $LOCAL_PATH
    $MPIRUN -np $NPROC $SEM_PATH/bin/xsmooth_sem_sph_pde $sigma_h $sigma_v $param $LOCAL_PATH ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}.ls/ $GPU_MODE >> $logfile
    \rm $LOCAL_PATH/*_$param.bin
    for i in `seq 1 $NPROC`;
    do
      ii=`echo $i |awk '{printf "%06d", $1-1}'`
      name=${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}.ls/proc${ii}_${param}
      mv ${name}_smooth.bin $name.bin 
    done

    echo "converting $param to hdf5 ..."
    fwat-main bin2h5 ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}.ls/ $param $NPROC 1
    \rm ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}.ls/*_${param}.bin
  done 

  for param in hess_kernel;
  do 
    echo "converting $param to hdf5 ..."
    fwat-main bin2h5 ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}.ls/ $param $NPROC 1
    \rm ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}.ls/*_${param}.bin
  done
  echo " " 

  # check wolfe condition
  echo "line search ..."
  $MPIRUN $hostfile -np $NPROC  fwat-main linesearch $MODEL $chi $chi1 

  # check if this line search is accepted
  LSDIR=./${FWAT_OPT_DIR}/MODEL_${MODEL}.ls
  flag=`fwat-utils getparam flag ${LBFGS_FILE}`
  if [ "$flag" == "GRAD" ]; then 
    icur=$(echo $MODEL |awk -F'M' '{print $2}')
    inext=$(printf "%02d" `echo $MODEL |awk -F'M' '{print $2+1}'`)
    echo misfit for iteration $icur and $inext $chi $chi1 >> misfit.log

    echo " " >> $logfile
    echo "rename  MODEL_${MODEL}.ls =>  MODEL_M$inext" >> $logfile
    
    # move new model to optimize/MODEL_M$inext
    rm -rf ./${FWAT_OPT_DIR}/MODEL_M$inext 
    mv $LSDIR ./${FWAT_OPT_DIR}/MODEL_M$inext 

    # move kernels to optimize/SUM_KERNELS_M$inext
    rm -rf ./${FWAT_OPT_DIR}/SUM_KERNELS_M$inext
    mv ./${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}.ls ./${FWAT_OPT_DIR}/SUM_KERNELS_M$inext
    
    # solver 
    rm -rf ${FWAT_SOLVER}/M$inext
    mv ${FWAT_SOLVER}/${MODEL}.ls ${FWAT_SOLVER}/M$inext

    # misfits
    rm -rf ${FWAT_MISFIT}/M$inext 
    mv ${FWAT_MISFIT}/${MODEL}.ls ${FWAT_MISFIT}/M$inext 

    # save LOGS
    mkdir -p LOG/$MODEL LOG/M$inext
    cd LOG
    for f in  ADJ* POST* output_fwat[1,2]*;
    do 
      if [   -f $f ]; then 
        mv $f $MODEL/
      fi
    done 
    for f in LS* WOLFE* output_fwat[3,4]*;
    do 
      if [   -f $f ]; then 
        mv $f M$inext/
      fi
    done 
    cd ..

    # clean useless information
    echo " " >> $logfile
    echo "clean useless files" >> $logfile
    for d in $MODEL M$inext;
    do 
      for CDIR in SEM GRADIENT;do 
        for f in ${FWAT_SOLVER}/$d/*/*$CDIR;
        do 
          echo "clean $f" >> $logfile 
          rm -rf $f 
        done 
      done
    done 

    echo " Finish line search direction  here!!!" >> $logfile 
  else 
    echo "$MPIRUN $hostfile -np $NPROC fwat-main update $MODEL $LSDIR"
    $MPIRUN $hostfile -np $NPROC fwat-main update $MODEL $LSDIR >> $logfile

    # generate new model database
    change_par LOCAL_PATH $LSDIR DATA/Par_file
    change_par LOCAL_PATH $LSDIR  DATA/meshfem3D_files/Mesh_Par_file
    change_par SAVE_MESH_FILES .false. DATA/Par_file

    # copy info to new 
    LOCAL_PATH=./${FWAT_OPT_DIR}/MODEL_${MODEL}
    echo -e ".false.\n.true." > adepml_stage
    \cp  $LOCAL_PATH/*Database $LSDIR/
    \cp  $LOCAL_PATH/*adepml* $LSDIR/
    \cp  $LOCAL_PATH/*undeformed_xyz.bin $LSDIR/
    :> DATA/FORCESOLUTION
    $MPIRUN $hostfile -np $NPROC $SEM_PATH/bin/xgenerate_databases 

    \rm adepml_*

    step_fac=`fwat-utils getparam alpha ${LBFGS_FILE}`
    echo " Line search failed, try step_fac = $step_fac !!!" >> $logfile 
  fi 

  echo " " >> $logfile

}

run_post () { # run post-processing
  # create log file
  logfile=LOG/output_fwat2_post_log_${MODEL}.txt
  :> $logfile
  echo "running POST " >> $logfile 

  # get smooth parameters
  GPU_MODE=`grep ^"GPU_MODE" DATA/Par_file | cut -d'=' -f2`
  LOCAL_PATH=${FWAT_OPT_DIR}/MODEL_${MODEL}
  change_par LOCAL_PATH $LOCAL_PATH ./DATA/Par_file
  change_par LOCAL_PATH $LOCAL_PATH ./DATA/meshfem3D_files/Mesh_Par_file
  info=`fwat-utils getparam optimize/SMOOTHING  | sed 's/\[\|]//g' | sed 's/,/ /g'`
  sigma_h=`echo $info | awk  '{print $1}'`
  sigma_v=`echo $info | awk  '{print $2}'`

  # sum kernels
  if [ $FLAG != "GRAD" ]; then 
    echo "sum kernels ..."
    echo "CMD: fwat-main sum_kernel  $MODEL"
    $MPIRUN $hostfile -np $NPROC fwat-main sum_kernel $MODEL >> $logfile
    echo " " >> $logfile

    # smooth gradient 
    kl_list=`fwat-model name direc`
    for param1 in $kl_list; 
    do 
      param=${param1:1}_kernel
      mv ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}/*_$param.bin $LOCAL_PATH
      $MPIRUN -np $NPROC $SEM_PATH/bin/xsmooth_sem_sph_pde $sigma_h $sigma_v $param $LOCAL_PATH ${FWAT_OPT_DIR}/SUM_KERNELS_$MODEL/ $GPU_MODE >> $logfile
      \rm $LOCAL_PATH/*_$param.bin
      for i in `seq 1 $NPROC`;
      do
        ii=`echo $i |awk '{printf "%06d", $1-1}'`
        name=${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}/proc${ii}_${param}
        mv ${name}_smooth.bin $name.bin 
      done

      echo "converting $param to hdf5 ..."
      fwat-main bin2h5 ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}/ $param $NPROC 1
      \rm ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}/*_${param}.bin
    done
  fi 

  echo " "

  # smooth hess kernel if required
  if [ $PRECOND == "default" ] && [ $MODEL == "M00"  ];then 
    param=hess_kernel
    mv ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}/*_$param.bin $LOCAL_PATH
    $MPIRUN -np $NPROC $SEM_PATH/bin/xsmooth_sem_sph_pde 50000 25000 $param $LOCAL_PATH ${FWAT_OPT_DIR}/SUM_KERNELS_$MODEL/ $GPU_MODE >> $logfile
    \rm $LOCAL_PATH/*_$param.bin
    for i in `seq 1 $NPROC`;
    do
      ii=`echo $i |awk '{printf "%06d", $1-1}'`
      name=${FWAT_OPT_DIR}/SUM_KERNELS_$MODEL/proc${ii}_$param
      mv ${name}_smooth.bin $name.bin 
    done 
    \rm ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}/*_${param}.bin
  fi 

  if [ $FLAG != "GRAD" ]; then 
    fwat-main bin2h5 ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}/ hess_kernel $NPROC 1
    \rm ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}/*hess_kernel.bin
  fi

  # get search direction
  $MPIRUN -np $NPROC fwat-main direc 
  echo " "

  # get search direction
  kl_list=`fwat-model name direc`
  for param in $kl_list; 
  do 
    echo "converting $param to hdf5 ..."
    fwat-main bin2h5 ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}/ $param $NPROC 1
    \rm ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}/*_${param}.bin
  done

  # generate new model
  LSDIR=./${FWAT_OPT_DIR}/MODEL_${MODEL}.ls
  mkdir -p $LSDIR
  echo " "
  echo "$MPIRUN -np $NPROC fwat-main update $MODEL $LSDIR"
  $MPIRUN -np $NPROC fwat-main update $MODEL $LSDIR >> $logfile

  # generate new model database
  change_par LOCAL_PATH $LSDIR DATA/Par_file
  change_par LOCAL_PATH $LSDIR  DATA/meshfem3D_files/Mesh_Par_file
  change_par SAVE_MESH_FILES .false. DATA/Par_file

  # copy info to new 
  echo -e ".false.\n.true." > adepml_stage
  \cp  $LOCAL_PATH/*Database $LSDIR/
  \cp  $LOCAL_PATH/*adepml* $LSDIR/
  \cp  $LOCAL_PATH/*undeformed_xyz.bin $LSDIR/
  :> DATA/FORCESOLUTION
  $MPIRUN -np $NPROC $SEM_PATH/bin/xgenerate_databases 

  # delete 
  \rm adepml_*

  # check if search method is GD
  OPT_METHOD=`fwat-utils getparam optimize/OPT_METHOD`
  if [ "$OPT_METHOD" == "GD" ]; then 
    let iter1=iter+1
    fwat-utils setparam flag INIT ${LBFGS_FILE}
    fwat-utils setparam iter_start $iter1 ${LBFGS_FILE}
    fwat-utils setparam iter $iter1 ${LBFGS_FILE}
    fwat-utils setparam alpha -1. ${LBFGS_FILE}
    MODEL1=M`echo "$iter" |awk '{printf "%02d",$1+1}'`
    mv $LSDIR ./${FWAT_OPT_DIR}/MODEL_${MODEL1}

    # save LOGS
    cd LOG
    mkdir -p $MODEL
    for f in  ADJ* POST* output_fwat[1,2]*;
    do 
      if [   -f $f ]; then 
        mv $f $MODEL/
      fi
    done
    cd ..

    # clean useless information
    echo " " >> $logfile
    echo "clean useless files" >> $logfile
    for d in $MODEL;
    do 
      for CDIR in SEM GRADIENT;do 
        for f in ${FWAT_SOLVER}/$d/*/*$CDIR;
        do 
          echo "clean $f" >> $logfile 
          rm -rf $f 
        done 
      done
    done 
  fi

  echo " " >> $logfile
  echo "******************************************************" >> $logfile
  echo " Finished FWAT POST here!!!" >> $logfile 
  echo " " >> $logfile
}

if [ $# -ne 1 ]; then
  echo "Usage: bash sbash_post.sh post/wolfe"
  exit 1
fi
POSTID=$1
if [ "$POSTID" != "post" ] && [ "$POSTID" != "wolfe" ]; then
  echo "Usage: bash sbash_post.sh post/wolfe"
  exit 1
fi

$hostfile

# load parameters 
source config.env
source utils.sh

# get search direction
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`
iter=`fwat-utils getparam iter ${LBFGS_FILE}`
FLAG=`fwat-utils getparam flag ${LBFGS_FILE}`
MODEL=M`echo "$iter" |awk '{printf "%02d",$1}'`
PRECOND=`fwat-utils getparam ${FWAT_OPT_DIR}/PRECOND_TYPE`
SIMU_TYPES=(`fwat-utils getparam simulation/types| tr -d '[]",'\'`)
SIMU_WEIGHTS=(`fwat-utils getparam simulation/weights| tr -d '[]",'\'`)

# create a hostfile to run mpi job if required
GLOBAL_SLOTS="all_slots.txt"
hostfile=""
if [  -f "$GLOBAL_SLOTS" ]; then
  if [ -n "$SLURM_JOB_ID" ] || [ -n "$PBS_JOBID" ]; then
    head -n $NPROC $GLOBAL_SLOTS > hostfile_$POSTID.txt
    hostfile="--hostfile hostfile_$POSTID.txt"
  fi
fi

# SIMU_TYPES is defined in parameters.sh, which is a list of simulation types to run, e.g., SIMU_TYPES=( "fullwaveform" "travel_time" )
nsimtypes="${#SIMU_TYPES[@]}" 

if [ "$POSTID" == "post" ]; then
  run_post
elif [ "$POSTID" == "wolfe" ]; then
  run_wolfe
else 
  echo "Usage: bash sbash_post.sh post/wolfe"
  exit 1
fi

if [  -f "$GLOBAL_SLOTS" ]; then
  # clean hostfile
  if [ -n "$SLURM_JOB_ID" ] || [ -n "$PBS_JOBID" ]; then
    \rm hostfile_$POSTID.txt
  fi
fi


