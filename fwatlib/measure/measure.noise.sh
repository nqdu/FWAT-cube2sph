#!/bin/bash
if [[ $# -ne 3 ]]; then 
    echo "Usage: ./measure_tele.sh iter evtid run_opt(1,2,3)"
    exit 1
fi

source module_env
. parameters.sh

set -e

# get input args
iter=$1
evtid=$2
run_opt=$3
lsflag=""
if  [ $run_opt -eq 2 ]; then 
  lsflag=".ls"
fi 

ncpu=`cat /proc/cpuinfo |grep cores |wc -l |awk '{print $1/2}'`

# set directory
current_dir=`pwd`
mod=M`printf "%02d" $iter`
rundir=solver/$mod"$lsflag"/$evtid/
echo "convert synthetic seismograms to sac : evtid=$evtid iter=$iter ..."
SYN_DIR=$rundir/OUTPUT_FILES 

# rotate seismograms to ENZ
mpirun -np 4 python  \
        $MEASURE_LIB/rotate_seismogram.py --fn_matrix="src_rec/rot_$evtid"   \
       --rotate="XYZ->NEZ" --from_dir="$rundir/OUTPUT_FILES/" --to_dir="$rundir/OUTPUT_FILES/"   \
       --from_template='${nt}.${sta}.BX${comp}.semd'  \
       --to_template='${nt}.${sta}.BX${comp}.sem.ascii'

# ascii to sac
echo " python $MEASURE_LIB/ascii2sac.py $current_dir/src_rec/STATIONS_${evtid}_globe  \
      $current_dir/src_rec/FORCE_ORG/FORCESOLUTION_${evtid} $SYN_DIR"
\rm -f $SYN_DIR/*.sac # clean all sac file
python $MEASURE_LIB/ascii2sac.py $current_dir/src_rec/STATIONS_${evtid}_globe  \
      $current_dir/src_rec/FORCE_ORG/FORCESOLUTION_${evtid} $SYN_DIR
\rm -f $SYN_DIR/*.sem.ascii

# run measure and write measure_adj file
echo " "
if [ $run_opt -ne 1 ]; then 
  mpirun -np $ncpu python $MEASURE_LIB/preprocess_noise.py $iter $evtid $run_opt
else 
  mpirun -np $ncpu python $MEASURE_LIB/preprocess_noise.py $iter $evtid $run_opt
  exit 0
fi 

# get filter params
# measure misifts
cd $SYN_DIR
fwat_file=../DATA/FWAT.PAR.yaml
info=`python $FWATLIB/get_param.py measure/noise/FILTER_BANDS  $fwat_file| sed 's/\[\|]//g' | sed 's/,/ /g'`
LONG_P=(`echo $info | awk '{for(i=2; i<=NF; i+=2) print $i}'`)
SHORT_P=(`echo $info | awk '{for(i=1; i<=NF; i+=2) print $i}'`)
NUM_FILTER=`echo ${#SHORT_P[@]}`
mkdir -p $current_dir/misfits/
outfile=$current_dir/output_fwat1_log.$mod.$evtid"$lsflag".txt 
if [ $run_opt == 2 ]; then
  outfile=$current_dir/output_fwat3_log.$mod.$evtid"$lsflag".txt 
fi 

# run measure_adj 
echo " "
for ((i=0;i<$NUM_FILTER;i++));
do
  band=`printf "T%03g_T%03g" ${SHORT_P[$i]} ${LONG_P[$i]}`
  echo "measure adj for $band ..."
  cd $band 
  mkdir -p OUTPUT_FILES
  \cp ../../DATA/MEASUREMENT.PAR  .
  mpirun -np $ncpu $MEASURE_LIB/bin/measure_adj_mpi
  awk -v a=$evtid '{$1=a;$29=$29*1.; print}' window_chi >  $current_dir/misfits/$mod.${evtid}${lsflag}_${band}_noise_window_chi
  
  # # continue for line search
  # if [ $run_opt == 2 ]; then
  #   \rm -rf OUTPUT_FILES
  #   cd $current_dir
  #   cd $SYN_DIR
  #   continue
  # fi

  # convert sac to hdf5
  python $FWATLIB/measure/sac2h5.py ../seismogram.obs.$band.h5 *.sac.obs
  python $FWATLIB/measure/sac2h5.py ../seismogram.syn.$band.h5 *.sac.syn

  # kerstr
  kerstr=iker`head -2  MEASUREMENT.PAR  |tail -1 | awk '{printf "%02d", $1}'`

  # rotate adjoint source 
  cd OUTPUT_FILES
  for f in `ls *iker*.adj |cut -d'.' -f1,2 |sort -n |uniq`;
  do 
    newf=`echo $f |awk -F'.' '{new_var=$2"."$1; print new_var}'`
    baz=`saclst baz f ../$newf.BXZ.sac.obs | awk '{print $2 * 3.1415926535898/180.}'`

    $MEASURE_LIB/bin/rotate_adj_src $baz $f.BXZ.$kerstr.adj $f.BXT.$kerstr.adj $f.BXR.$kerstr.adj $newf.BXE.adj $newf.BXN.adj 
    \cp $f.BXZ.$kerstr.adj $newf.BXZ.adj
  done
  cd ../..

done

# # continue for finish line search
# if [ $run_opt -eq 2 ]; then
#   exit 0;
# fi

cd $current_dir
bash $MEASURE_LIB/sum_adj_source.sh $iter $evtid noise $run_opt
