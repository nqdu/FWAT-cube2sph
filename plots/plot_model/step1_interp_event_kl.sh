#!/bin/bash
source parameters.sh
source module_env

set -e 

# work for ls model

########## stop ##################

param_set=$KERNEL_SET
WORK_PREFIX=$SOLVER_DIR/M

lsflag=""
if [ "$INTP_LS" == "1" ]; then 
  lsflag=".ls"
fi 

# interpolate
echo ""
echo "interpolating ..."
for iter in $run_indx;do 
  ii=`printf %02d $iter`
  # get how many events
  for evt in `ls ${WORK_PREFIX}${ii}${lsflag}/`;
  do 
    workdir=${WORK_PREFIX}${ii}${lsflag}/$evt/GRADIENT
    $MPIRUN -np 4 python $fwatlib/optimize/write_event_kernels.py $MODEL_DIR/MODEL_M${ii}${lsflag} $workdir $MDTYPE $KLTYPE

    for param in $param_set ;do 
    for name in horiz verti ;  do
      nfiles=`ls input/ |grep $name.*.loc |wc -l`
      for ip in `seq 1 $nfiles`; do
        # interpolate
        $specfem_dir/bin/xcreate_slice $param $workdir  $DATABASE_DIR \
            input/$name.$ip.txt input/$name.$ip.loc $param.$evt.$name.$ip.$ii.out .false.
      done
    done
    done

    # remove binaries
    \rm $workdir/*.bin
  done
done

echo ""
echo "create txt for ploting ..."
mkdir -p profiles/

for param in $param_set ;do 
for iter in $run_indx;do 
  ii=`printf %02d $iter`
  for evt in `ls ${WORK_PREFIX}${ii}${lsflag}/`;
  do
    name=horiz
    nfiles=`ls input/ |grep $name.*.loc |wc -l`
    for ip in `seq 1 $nfiles`; do
    
      
      awk '{print $4,$5}' input/$name.$ip.txt > tmp.1 
      awk '{print $1}' $param.$evt.$name.$ip.$ii.out > tmp.2 
      paste tmp.1 tmp.2 > tmp.3 
      mv tmp.3 profiles/$param.$evt.iter$ii${lsflag}.$name.$ip.txt
      \rm tmp.1 tmp.2 $param.$evt.$name.$ip.$ii.out
    done

    name=verti
    nfiles=`ls input/ |grep $name.*.loc |wc -l`
    for ip in `seq 1 $nfiles`; do
      awk '{print $5,-$4/1000}' input/$name.$ip.txt > tmp.1 
      awk '{print $1}' $param.$evt.$name.$ip.$ii.out > tmp.2 
      paste tmp.1 tmp.2 > tmp.3 
      mv tmp.3 profiles/$param.$evt.iter${ii}${lsflag}.$name.$ip.txt
      \rm tmp.1 tmp.2 $param.$evt.$name.$ip.$ii.out
    done
  done
done
done

# interpolate event kernels 