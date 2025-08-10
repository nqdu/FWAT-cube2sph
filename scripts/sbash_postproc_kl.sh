set -e 
# include file
source module_env
source parameters.sh

# input vars
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`

# get search direction
FWATPARAM=./fwat_params
iter=`fwat-utils getparam iter $FWATPARAM/lbfgs.yaml`
FLAG=`fwat-utils getparam flag $FWATPARAM/lbfgs.yaml`
MODEL=M`echo "$iter" |awk '{printf "%02d",$1}'`
PRECOND=`fwat-utils getparam optimize/PRECOND_TYPE`

# check how many simu types required
nsimtypes="${#SIMU_TYPES[@]}"
if [ "$nsimtypes" == "1" ]; then 
  SOURCE_FILE=./src_rec/sources.dat.${SIMU_TYPES[0]}
else
  # generate files if requireds
  iter_start=`fwat-utils getparam iter_start $FWATPARAM/lbfgs.yaml`
  if [ "$iter" == "0" ]; then
    # init source file
    SOURCE_FILE=./src_rec/sources.dat.joint
    cat ./src_rec/sources.dat.${SIMU_TYPES[0]} > $SOURCE_FILE

    # compute misfit 
    info=`fwat-main misfit $MODEL ${SIMU_TYPES[0]}`
    chi0=`echo $info |awk '{print $1/$2}'`
    chi1=`echo $chi0 $chi0 ${SIMU_TYPES_USER_WEIGHT[0]} |awk '{print $1/$2*$3}'`

    # init weight_kl.txt 
    :> ./optimize/weight_kl.txt
    awk -v a=$chi1 '{print $2*0+a}' ./src_rec/sources.dat.${SIMU_TYPES[0]} >> ./optimize/weight_kl.txt

    # for other simulation types
    # misfit = L0 + L1 * s0/s1 + L2 * s0/s2 + L_i s0/s_i
    for((i=1;i<$nsimtypes;i++)); 
    do 
      cat ./src/sources.dat.${SIMU_TYPES[$i]} >> $SOURCE_FILE
      info=`fwat-main misfit $MODEL ${SIMU_TYPES[$i]}`
      chi1=`echo $info |awk '{print $1/$2}'`
      chi1=`echo $chi0 $chi1 ${SIMU_TYPES_USER_WEIGHT[$i]} |awk '{print $1/$2*$3}'`
      awk -v a=$chi1 '{print $2*0+a}' ./src_rec/sources.dat.${SIMU_TYPES[$i]} >> ./optimize/weight_kl.txt
    done 
  fi
fi 

# create log file
logfile=LOG/output_fwat2_post_log_${MODEL}.txt
:> $logfile
echo "running POST " >> $logfile 

# sum kernels
if [ $FLAG != "GRAD" ]; then 
  echo "sum kernels ..."
  echo "CMD: fwat-main sum_kernel $SOURCE_FILE $iter $MODEL"
  $MPIRUN -np $NPROC fwat-main sum_kernel $SOURCE_FILE $iter $MODEL >> $logfile

  kl_list=`fwat-model name grad`
  for param in $kl_list 
  do 
    echo "converting $param to hdf5 ..."
    fwat-main bin2h5 optimize/SUM_KERNELS_${MODEL}/ $param $NPROC 1
    \rm optimize/SUM_KERNELS_${MODEL}/*_${param}.bin
  done
  echo " " 
fi 

# set smoothing parameters
LOCAL_PATH=./optimize/MODEL_${MODEL}
change_par LOCAL_PATH $LOCAL_PATH ./DATA/Par_file
change_par LOCAL_PATH $LOCAL_PATH ./DATA/meshfem3D_files/Mesh_Par_file
info=`fwat-utils getparam optimize/SMOOTHING  | sed 's/\[\|]//g' | sed 's/,/ /g'`
sigma_h=`echo $info | awk  '{print $1}'`
sigma_v=`echo $info | awk  '{print $2}'`

# smooth hess kernel if required
GPU_MODE=`grep ^"GPU_MODE" DATA/Par_file | cut -d'=' -f2`
if [ $PRECOND == "default" ] && [ $MODEL == "M00"  ];then 
  param=hess_kernel
  mv optimize/SUM_KERNELS_${MODEL}/*_$param.bin $LOCAL_PATH
  $MPIRUN -np $NPROC $SEM_PATH/bin/xsmooth_sem_sph_pde 50000 25000 $param $LOCAL_PATH optimize/SUM_KERNELS_$MODEL/ $GPU_MODE >> $logfile
  \rm $LOCAL_PATH/*_$param.bin
  for i in `seq 1 $NPROC`;
  do
    ii=`echo $i |awk '{printf "%06d", $1-1}'`
    name=optimize/SUM_KERNELS_$MODEL/proc${ii}_$param
    mv ${name}_smooth.bin $name.bin 
  done 
  \rm optimize/SUM_KERNELS_${MODEL}/*_${param}.bin
fi 

if [ $FLAG != "GRAD" ]; then 
  fwat-main bin2h5 optimize/SUM_KERNELS_${MODEL}/ hess_kernel $NPROC 1
  \rm optimize/SUM_KERNELS_${MODEL}/*hess_kernel.bin
fi

# get search direction
$MPIRUN -np $NPROC fwat-main direc 
echo " "

# smooth search direction
kl_list=`fwat-model name direc`
for param in $kl_list; 
do 
  mv optimize/SUM_KERNELS_$MODEL/*_$param.bin $LOCAL_PATH
  $MPIRUN -np $NPROC $SEM_PATH/bin/xsmooth_sem_sph_pde $sigma_h $sigma_v $param $LOCAL_PATH optimize/SUM_KERNELS_$MODEL/ $GPU_MODE >> $logfile
  \rm $LOCAL_PATH/*_$param.bin
  for i in `seq 1 $NPROC`;
  do
    ii=`echo $i |awk '{printf "%06d", $1-1}'`
    name=optimize/SUM_KERNELS_${MODEL}/proc${ii}_$param
    mv ${name}_smooth.bin $name.bin 
  done

  echo "converting $param to hdf5 ..."
  fwat-main bin2h5 optimize/SUM_KERNELS_${MODEL}/ $param $NPROC 1
  \rm optimize/SUM_KERNELS_${MODEL}/*_${param}.bin
done

# generate new model
LSDIR=./optimize/MODEL_${MODEL}.ls
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
$MPIRUN -np $NPROC $SEM_PATH/bin/xgenerate_databases 

# delete 
\rm adepml_*

# check if search method is GD
OPT_METHOD=`fwat-utils getparam optimize/OPT_METHOD`
if [ "$OPT_METHOD" == "GD" ]; then 
	fwat-utils setparam flag INIT $FWATPARAM/lbfgs.yaml
	fwat-utils setparam iter_start $iter $FWATPARAM/lbfgs.yaml
	let iter1=iter+1
	fwat-utils setparam iter $iter1 $FWATPARAM/lbfgs.yaml
	MODEL1=M`echo "$iter" |awk '{printf "%02d",$1+1}'`
	mv $LSDIR ./optimize/MODEL_${MODEL1}

	# save LOGS
	cd LOG
	mkdir LOG/$MODEL
  for f in  FWD_ADJ* POST* output_fwat[1,2]*;
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
      for f in solver/$d/*/*$CDIR;
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
