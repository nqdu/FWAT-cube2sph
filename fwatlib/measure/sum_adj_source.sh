
if [[ $# -ne 2 ]]; then 
    echo "Usage: ./sum_adj_source iter evtid"
    exit 1
fi

source module_env
. parameters.sh


set -e

# get input args
iter=$1
evtid=$2
lsflag=""

# set directory
current_dir=`pwd`
mod=M`printf "%02d" $iter`
rundir=solver/$mod"$lsflag"/$evtid/
SYN_DIR=$rundir/OUTPUT_FILES 

# get frequency band
cd $SYN_DIR
fwat_file=../DATA/FWAT.PAR
SHORT_P=(`cat $fwat_file |grep 'SHORT_P:' |awk -F: '{print $2}'`)
LONG_P=(`cat $fwat_file |grep 'LONG_P:' |awk -F: '{print $2}'`)
NUM_FILTER=`echo ${#SHORT_P[@]}`

echo " "
echo "computing ascii adjoint source ..." 
cd $current_dir
rm -rf $rundir/SEM; mkdir -p $rundir/SEM
cd $rundir/OUTPUT_FILES
band0=`printf "T%03d_T%03d" ${SHORT_P[0]} ${LONG_P[0]}`
filenames=`ls ${band0}/OUTPUT_FILES/ |grep iker | cut -d'.' -f1,2 |sort -n |uniq`
for f in $filenames
do 
  newf=`echo $f |awk -F'.' '{new_var=$2"."$1; print new_var}'`
  for c in Z N E;
  do
    for ((i=0;i<$NUM_FILTER;i++));
    do
      band=`printf "T%03d_T%03d" ${SHORT_P[$i]} ${LONG_P[$i]}`
      if [ $i -eq 0 ]; then 
        awk '{print $1,0.}' ${band0}/OUTPUT_FILES/$newf.BX$c.adj  > temp0
      fi 
      paste temp0 ${band}/OUTPUT_FILES/$newf.BX$c.adj |awk '{print $1,$2+$4}' > temp1 
      cat temp1 > temp0 
      rm temp1
    done

    mv temp0 ../SEM/$newf.BX$c.adj.ascii
  done 
done

cd $current_dir
# rotate seismograms to xyz
echo "--fn_matrix="$current_dir/src_rec/rot_$evtid"   \
       --rotate='XYZ<-NEZ' --from_dir="$rundir/SEM" --to_dir="$rundir/SEM""
mpirun -np 4 python $MEASURE_LIB/rotate_seismogram.py --fn_matrix="$current_dir/src_rec/rot_$evtid"   \
       --rotate='XYZ<-NEZ' --from_dir="$rundir/SEM" --to_dir="$rundir/SEM"   \
       --from_template='${nt}.${sta}.BX${comp}.adj.ascii'  \
       --to_template='${nt}.${sta}.BX${comp}.adj'
echo ""
