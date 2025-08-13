#!/bin/bash
source parameters.sh
source module_env

# work for ls model

########## stop ##################

if [ "$INTP_KL" == "1" ]; then 
  param_set=$KERNEL_SET
  WORK_PREFIX=$MODEL_DIR/SUM_KERNELS_M
else 
  param_set=$MODEL_SET
  WORK_PREFIX=$MODEL_DIR/MODEL_M
fi

lsflag=""
if [ "$INTP_LS" == "1" ]; then 
  lsflag=".ls"
fi 

# interpolate
echo ""
echo "interpolating ..."
for param in $param_set ;do 
for iter in $run_indx;do 
  ii=`printf %02d $iter`
  for name in horiz verti ;  do
    nfiles=`ls input/ |grep $name.*.loc |wc -l`
    for ip in `seq 1 $nfiles`; do

      if [ "$INTP_KL" == "1" ]; then 
        fwat-main bin2h5 ${WORK_PREFIX}${ii}${lsflag}/  $param $NPROC 0
      fi

      $specfem_dir/bin/xcreate_slice $param ${WORK_PREFIX}${ii}${lsflag}/  $DATABASE_DIR \
          input/$name.$ip.txt input/$name.$ip.loc $param.$name.$ip.$ii.out .false.
      
      if [ "$INTP_KL" == "1" ]; then 
        \rm ${WORK_PREFIX}$ii/*$param.bin
      fi
    done
  done
done
done

echo ""
echo "create txt for ploting ..."
mkdir -p profiles/

for param in $param_set ;do 
for iter in $run_indx;do 
  name=horiz
  for ip in `seq 1 $NSLICE_HORIZ`; do
  
    ii=`printf %02d $iter`
    awk '{print $4,$5}' input/$name.$ip.txt > tmp.1 
    awk '{print $1}' $param.$name.$ip.$ii.out > tmp.2 
    paste tmp.1 tmp.2 > tmp.3 
    mv tmp.3 profiles/$param.iter$ii${lsflag}.$name.$ip.txt
    \rm tmp.1 tmp.2 $param.$name.$ip.$ii.out
  done

  name=verti
  for ip in `seq 1 $NSLICE_VERTI`; do
    ii=`printf %02d $iter`
    awk '{print $5,-$4/1000}' input/$name.$ip.txt > tmp.1 
    awk '{print $1}' $param.$name.$ip.$ii.out > tmp.2 
    paste tmp.1 tmp.2 > tmp.3 
    mv tmp.3 profiles/$param.iter${ii}${lsflag}.$name.$ip.txt
    \rm tmp.1 tmp.2 $param.$name.$ip.$ii.out
  done
done
done