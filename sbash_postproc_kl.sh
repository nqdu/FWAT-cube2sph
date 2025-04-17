#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=160
#SBATCH --time=00:20:59
#SBATCH --job-name=POST
#SBATCH --output=POST_%j.txt
#SBATCH --partition=compute 

set -e 
# include file
. parameters.sh
# load modules 
source module_env

# input vars
SOURCE_FILE=./src_rec/sources.dat.tele
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`

# parfile changer script
change_par=$FWATLIB/change_par_file.sh

# get search direction
iter=`python $FWATLIB/get_param.py iter $FWATPARAM/lbfgs.yaml`
FLAG=`python $FWATLIB/get_param.py flag $FWATPARAM/lbfgs.yaml`
MODEL=M`echo "$iter" |awk '{printf "%02d",$1}'`
PRECOND=`python $FWATLIB/get_param.py optimize/PRECOND_TYPE`


# create log file
logfile=output_fwat2_post_log_${MODEL}.txt
:> $logfile
echo "running POST " >> $logfile 

# sum kernels
if [ $FLAG != "GRAD" ]; then 
  echo "sum kernels ..."
  echo "CMD: mpirun -np $NPROC python $OPT_LIB/sum_kernel.py $SOURCE_FILE $iter $PRECOND $MODEL"
  mpirun -np $NPROC python $OPT_LIB/sum_kernel.py $SOURCE_FILE $iter $PRECOND $MODEL >> $logfile

  kl_list=`GET_GRAD_NAME`
  for param in $kl_list 
  do 
    echo "converting $param to hdf5 ..."
    python $OPT_LIB/bin2h5.py optimize/SUM_KERNELS_${MODEL}/ $param $NPROC 1
    \rm optimize/SUM_KERNELS_${MODEL}/*_${param}.bin
  done
  echo " " 
fi 

# set smoothing parameters
LOCAL_PATH=./DATABASES_MPI
$change_par LOCAL_PATH $LOCAL_PATH ./DATA/Par_file
$change_par LOCAL_PATH $LOCAL_PATH ./DATA/meshfem3D_files/Mesh_Par_file
info=`python $FWATLIB/get_param.py optimize/SMOOTHING  | sed 's/\[\|]//g' | sed 's/,/ /g'`
sigma_h=`echo $info | awk  '{print $1}'`
sigma_v=`echo $info | awk  '{print $2}'`

# smooth hess kernel if required
if [ $PRECOND == "default" ] && [ $MODEL == "M00"  ];then 
  param=hess_kernel
  mv optimize/SUM_KERNELS_${MODEL}/*_$param.bin $LOCAL_PATH
  mpirun -np $NPROC $fksem/bin/xsmooth_sem_sph_pde 50000 25000 $param $LOCAL_PATH optimize/SUM_KERNELS_$MODEL/ .false.
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
  python $OPT_LIB/bin2h5.py optimize/SUM_KERNELS_${MODEL}/ hess_kernel $NPROC 1
  \rm optimize/SUM_KERNELS_${MODEL}/*hess_kernel.bin
fi

# get search direction
mpirun -np $NPROC python $OPT_LIB/get_lbfgs_direc.py $iter $FWATPARAM/lbfgs.yaml 
echo " "

# smooth search direction
kl_list=`GET_DIREC_NAME`
for param in $kl_list; 
do 
  mv optimize/SUM_KERNELS_$MODEL/*_$param.bin $LOCAL_PATH
  mpirun -np $NPROC $fksem/bin/xsmooth_sem_sph_pde $sigma_h $sigma_v $param $LOCAL_PATH optimize/SUM_KERNELS_$MODEL/ .false.
  \rm $LOCAL_PATH/*_$param.bin
  for i in `seq 1 $NPROC`;
  do
    ii=`echo $i |awk '{printf "%06d", $1-1}'`
    name=optimize/SUM_KERNELS_${MODEL}/proc${ii}_$param
    mv ${name}_smooth.bin $name.bin 
  done

  echo "converting $param to hdf5 ..."
  python $OPT_LIB/bin2h5.py optimize/SUM_KERNELS_${MODEL}/ $param $NPROC 1
  \rm optimize/SUM_KERNELS_${MODEL}/*_${param}.bin
done

# generate new model
LSDIR=./optimize/MODEL_${MODEL}.ls
mkdir -p $LSDIR
echo " "
echo "python $OPT_LIB/get_lbfgs_next.py $MODEL $LSDIR $FWATPARAM/FWAT.PAR.yaml $FWATPARAM/lbfgs.yaml  $NPROC"
python $OPT_LIB/get_lbfgs_next.py $MODEL $LSDIR $FWATPARAM/FWAT.PAR.yaml $FWATPARAM/lbfgs.yaml  $NPROC  >> $logfile

# generate new model database
$change_par LOCAL_PATH $LSDIR DATA/Par_file
$change_par LOCAL_PATH $LSDIR  DATA/meshfem3D_files/Mesh_Par_file
$change_par SAVE_MESH_FILES .false. DATA/Par_file

# copy info to new 
echo -e ".false.\n.true." > adepml_stage
\cp  $LOCAL_PATH/*Database $LSDIR/
\cp  $LOCAL_PATH/*adepml* $LSDIR/
\cp  $LOCAL_PATH/*undeformed_xyz.bin $LSDIR/
mpirun -np $NPROC $fksem/bin/xgenerate_databases 

# delete 
\rm adepml_*

echo " " >> $logfile
echo "******************************************************" >> $logfile
echo " Finished FWAT POST here!!!" >> $logfile 
echo " " >> $logfile
