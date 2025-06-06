
if [[ $# -ne 3 && $# -ne 4 ]]; then
    echo "Usage: ./sum_adj_source iter evtid simu_type (run_opt=3)"
    exit 1
fi

source module_env
. parameters.sh


set -e

# get input args
iter=$1
evtid=$2
simu_type=$3
run_opt=3
lsflag="" 
if [[ $# -eq 4 ]]; then
  run_opt=$4
  if [ $run_opt -eq 2 ]; then 
    lsflag=".ls"
  fi
fi

# set directory
current_dir=`pwd`
mod=M`printf "%02d" $iter`
rundir=solver/$mod"$lsflag"/$evtid/
SYN_DIR=$rundir/OUTPUT_FILES 

# get frequency band
cd $SYN_DIR
fwat_file=../DATA/FWAT.PAR.yaml
info=`python $FWATLIB/get_param.py measure/$simu_type/FILTER_BANDS $fwat_file | sed 's/\[\|]//g' | sed 's/,/ /g'`
LONG_P=(`echo $info | awk '{for(i=2; i<=NF; i+=2) print $i}'`)
SHORT_P=(`echo $info | awk '{for(i=1; i<=NF; i+=2) print $i}'`)
NUM_FILTER=`echo ${#SHORT_P[@]}`

echo " "
echo "computing ascii adjoint source ..." 
cd $current_dir
rm -rf $rundir/SEM; mkdir -p $rundir/SEM
cd $rundir/OUTPUT_FILES
band0=`printf "T%03g_T%03g" ${SHORT_P[0]} ${LONG_P[0]}`
filenames=`ls ${band0}/OUTPUT_FILES/ |grep adj | cut -d'.' -f1,2 |sort -n |uniq`
for f in $filenames
do 
  newf=`echo $f |awk -F'.' '{new_var=$1"."$2; print new_var}'`
  for c in Z N E;
  do
    paste T*_T*/OUTPUT_FILES/$newf.BX$c.adj | awk '{
  sum = 0;
  for (i=2; i<=NF; i+=2) sum += $i;
  print $1, sum
}' > ../SEM/$newf.BX$c.adj.ascii

  done 
done

cd $current_dir
# rotate seismograms to xyz
python $MEASURE_LIB/pack_seismogram.py $rundir/SEM/seismograms.h5 $rundir/SEM/*.adj.ascii
\rm $rundir/SEM/*.ascii
echo "--fn_matrix="$current_dir/src_rec/rot_$evtid"   \
       --rotate='XYZ<-NEZ' --from_dir="$rundir/SEM" --to_dir="$rundir/SEM""
mpirun -np $NPROC_PRE python $MEASURE_LIB/rotate_seismogram.py --fn_matrix="$current_dir/src_rec/rot_$evtid"   \
       --rotate='XYZ<-NEZ' --from_dir="$rundir/SEM" --to_dir="$rundir/SEM"   \
       --from_template='${nt}.${sta}.BX${comp}.adj.ascii'  \
       --to_template='${nt}.${sta}.BX${comp}.adj'
echo ""
