#!/bin/bash

if [[ $# -ne 3 ]]; then 
    echo "Usage: ./measure.sks.sh iter evtid run_opt(1,2,3)"
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

# set directory
current_dir=`pwd`
mod=M`printf "%02d" $iter`
rundir=solver/$mod"$lsflag"/$evtid/
echo "convert synthetic seismograms to sac : evtid=$evtid iter=$iter"
SYN_DIR=$rundir/OUTPUT_FILES 

# rotate sesmograms to ENZ
mpirun -np $NPROC_PRE python $MEASURE_LIB/rotate_seismogram.py --fn_matrix="src_rec/rot_$evtid"   \
       --rotate="XYZ->NEZ" --from_dir="$rundir/OUTPUT_FILES/" --to_dir="$rundir/OUTPUT_FILES/"   \
       --from_template='${nt}.${sta}.BX${comp}.semd'  \
       --to_template='${nt}.${sta}.BX${comp}.sem.ascii'

echo "python $MEASURE_LIB/ascii2sac.py $current_dir/src_rec/STATIONS_${evtid}_globe  \
      $current_dir/src_rec/CMTSOLUTION_${evtid} $SYN_DIR"
\rm -f $SYN_DIR/*.sac # clean all sac file
mpirun -np $NPROC_PRE python $MEASURE_LIB/ascii2sac.py $current_dir/src_rec/STATIONS_${evtid}_globe  \
      $current_dir/src_rec/CMTSOLUTION_${evtid} $SYN_DIR
\rm -f $SYN_DIR/*.sem.ascii

# rotate to RTZ
echo ""
echo "rotating ENZ TO RTZ ..."
cd $SYN_DIR
for name in `ls |grep E.sac |cut -d'.' -f1,2 |sort -n |uniq`;
do
  sac << EOF
r $name.BX[EN].sac 
rotate to gcp 
ch FILE 1 kcmpnm BHR
ch FILE 2 kcmpnm BHT
w $name.BXR.sac $name.BXT.sac
q
EOF
done 
cd $current_dir 

# preprocessing, miisfit, adjoint source 
echo " "
if [ $run_opt -ne 1 ]; then
  echo "estimating STF and synthetic data ..."
  mpirun -np $NPROC_PRE python $MEASURE_LIB/preprocess_sks.py $iter $evtid $run_opt
else
  echo "synthetic data ..."
  mpirun -np $NPROC_PRE python $MEASURE_LIB/preprocess_sks.py $iter $evtid $run_opt
  exit 0
fi

# measure misifts
echo ""
cd $SYN_DIR
fwat_file=../DATA/FWAT.PAR.yaml
info=`python $FWATLIB/get_param.py measure/tele/FILTER_BANDS $fwat_file  | sed 's/\[\|]//g' | sed 's/,/ /g'`
LONG_P=(`echo $info | awk '{for(i=2; i<=NF; i+=2) print $i}'`)
SHORT_P=(`echo $info | awk '{for(i=1; i<=NF; i+=2) print $i}'`)
NUM_FILTER=`echo ${#SHORT_P[@]}`
mkdir -p $current_dir/misfits/
outfile=$current_dir/output_fwat1_log.$mod.$evtid"$lsflag".txt 
if [ $run_opt == 2 ]; then
  outfile=$current_dir/output_fwat3_log.$mod.$evtid"$lsflag".txt 
fi 

for ((i=0;i<$NUM_FILTER;i++));
do
  band=`printf "T%03g_T%03g" ${SHORT_P[$i]} ${LONG_P[$i]}`
  cd $band 
  cat window_chi >  $current_dir/misfits/$mod.${evtid}${lsflag}_${band}_tele_window_chi

  # convert sac to hdf5
  python $FWATLIB/measure/sac2h5.py ../seismogram.obs.$band.h5 *.sac.obs
  python $FWATLIB/measure/sac2h5.py ../seismogram.syn.$band.h5 *.sac.syn
done

cd $current_dir
bash $MEASURE_LIB/sum_adj_source.sh $iter $evtid tele $run_opt