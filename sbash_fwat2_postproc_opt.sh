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
MODEL=M24
SOURCE_FILE=./src_rec/sources.dat.tele
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`

echo "running fwat2 " > output_fwat2_log_${MODEL}.txt

# parfile changer script
change_par=$FWATLIB/change_par_file.sh

# get search direction
iter=`echo "$MODEL" |awk -F'M' '{print $2}' | awk '{printf "%d",$1}'`

# copy model to optimize
LOCAL_PATH=./DATABASES_MPI
if [ $MODEL == "M00"  ];then
	rm -rf optimize/MODEL_M00; mkdir -p optimize/MODEL_M00
	\cp -r $LOCAL_PATH/* optimize/MODEL_M00
fi

# create log file
filename=output_fwat2_log_${MODEL}.txt
:> $filename

# sum kernels
echo "sum kernels ..."
PRECOND=`python $FWATLIB/get_param.py optimize/PRECOND_TYPE`
mpirun -np $NPROC python $OPT_LIB/sum_kernel.py  $iter $SOURCE_FILE $PRECOND >> $filename
kl_list=`GET_GRAD_NAME`
for param in $kl_list 
do 
  echo "\nconverting $param to hdf5 ..."
  python $OPT_LIB/bin2h5.py optimize/SUM_KERNELS_${MODEL}/ $param $NPROC 1
  \rm optimize/SUM_KERNELS_${MODEL}/*_${param}.bin
done
echo " " 

# set smoothing parameters
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
python $OPT_LIB/bin2h5.py optimize/SUM_KERNELS_${MODEL}/ hess_kernel $NPROC 1
\rm optimize/SUM_KERNELS_${MODEL}/*hess_kernel.bin

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

  echo "\nconverting $param to hdf5 ..."
  python $OPT_LIB/bin2h5.py optimize/SUM_KERNELS_${MODEL}/ $param $NPROC 1
  \rm optimize/SUM_KERNELS_${MODEL}/*_${param}.bin
done

# generate new model
LSDIR=./optimize/MODEL_${MODEL}_step01
mkdir -p $LSDIR
step_fac=`python $FWATLIB/get_param.py STEP_FAC fwat_params/lbfgs.yaml`
MIS_FILE=misfits/$MODEL.mis
echo " "
echo "python $OPT_LIB/get_lbfgs_step_fac.py $MODEL $LSDIR $step_fac $NPROC"
python $OPT_LIB/get_lbfgs_step_fac.py $MODEL $LSDIR $step_fac $NPROC  > $MIS_FILE
step_fac=`tail -1 $MIS_FILE |cut -d'=' -f2 |awk '{print $1}'`
dmax=`tail -1 $MIS_FILE |cut -d'=' -f2 |awk '{print $2}'`

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
echo " " >> $filename
echo "******************************************************" >> $filename
echo " Finished FWAT stage2 here!!!" >> $filename 
echo " " >> $filename
