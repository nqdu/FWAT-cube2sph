set -e 

# load parameters 
. parameters.sh 
source module_env
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`

# get search direction
iter=`fwat-utils getparam iter ${LBFGS_FILE}`
FLAG=`fwat-utils getparam flag ${LBFGS_FILE}`
MODEL=M`echo "$iter" |awk '{printf "%02d",$1}'`

# compute misfit
# check how many simu types required
nsimtypes="${#SIMU_TYPES[@]}"
if [ "$nsimtypes" == "1" ]; then 
  SOURCE_FILE_LS=./${FWAT_SRC_REC}/sources.dat.${SIMU_TYPES[0]}

  # compute misfits
  info=`fwat-main misfit $MODEL ${SIMU_TYPES[0]}`
  chi=`echo $info |awk '{print $1}'`
  info=`fwat-main misfit $MODEL.ls ${SIMU_TYPES[0]}`
  chi1=`echo $info |awk '{print $1}'`
else

  SOURCE_FILE_LS=./${FWAT_SRC_REC}/sources.dat.joint
  #init misifits
  chi=0.
  chi1=0.

  # get first model
  iter_start=`fwat-utils getparam iter_start ${LBFGS_FILE}`
  MSTART=M`echo "$iter_start" |awk '{printf "%02d",$1}'`
  MSTART=M00

  # for all simulation types, compute weighted misfits
  for((i=0;i<$nsimtypes;i++)); 
  do 
    info=`fwat-main misfit $MSTART ${SIMU_TYPES[$i]}`
    l0=`echo $info |awk '{print $1}'`
    info=`fwat-main misfit $MODEL ${SIMU_TYPES[$i]}`
    l1=`echo $info |awk '{print $1}'`
    info=`fwat-main misfit $MODEL.ls ${SIMU_TYPES[$i]}`
    l2=`echo $info |awk '{print $1}'`

    # weighted sum for iteration n
    # \sum_i L^(n)_i / L^(0)_i * L^(0)_0 * user_weight
    chi=`echo $chi $l1 $l0 ${SIMU_TYPES_USER_WEIGHT[$i]} |awk '{print $1+$2/$3*$4}'`
    chi1=`echo $chi1 $l2 $l0 ${SIMU_TYPES_USER_WEIGHT[$i]}|awk '{print $1+$2/$3*$4}'`
  done 

  # misfit for first model of dataset 1
  info=`fwat-main misfit $MSTART ${SIMU_TYPES[0]}`
  l0=`echo $info |awk '{print $1}'`
  chi=`echo $chi $l0 |awk '{print $1*$2}'`
  chi1=`echo $chi1 $l0 |awk '{print $1*$2}'`
fi 

echo "misfit current/next = $chi $chi1"
echo " "

# sum kernels for line search, save to optimize/sum_kernels_$MODEL.ls
echo "sum kernels for new model ..."
$MPIRUN -np $NPROC fwat-main sum_kernel $SOURCE_FILE_LS $iter $PRECOND $MODEL.ls
kl_list=`fwat-model name grad`
for param in $kl_list hess_kernel;
do 
    echo "converting $param to hdf5 ..."
    fwat-main bin2h5 ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}.ls/ $param $NPROC 1
    \rm ${FWAT_OPT_DIR}/SUM_KERNELS_${MODEL}.ls/*_${param}.bin
done
echo " " 

# check wolfe condition
echo "line search ..."
$MPIRUN -np $NPROC fwat-main linesearch $MODEL $chi $chi1 

logfile=LOG/output_fwat4_log_$MODEL.txt
echo "******************************************************" > $logfile

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
  echo "$MPIRUN -np $NPROC fwat-main update $MODEL $LSDIR"
  $MPIRUN -np $NPROC fwat-main update $MODEL $LSDIR >> $logfile

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
  $MPIRUN -np $NPROC $SEM_PATH/bin/xgenerate_databases 

  \rm adepml_*

  step_fac=`fwat-utils getparam alpha ${LBFGS_FILE}`
  echo " Line search failed, try step_fac = $step_fac !!!" >> $logfile 
fi 

echo " " >> $logfile
