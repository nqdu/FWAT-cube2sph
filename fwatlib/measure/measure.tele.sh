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

# set directory
current_dir=`pwd`
mod=M`printf "%02d" $iter`
rundir=solver/$mod"$lsflag"/$evtid/
echo "convert synthetic seismograms to sac : evtid=$evtid iter=$iter"
SYN_DIR=$rundir/OUTPUT_FILES 

# rotate seismograms to ENZ
mpirun -np 4 python $MEASURE_LIB/rotate_seismogram.py --fn_matrix="src_rec/rot_$evtid"   \
       --rotate="XYZ->NEZ" --from_dir="$rundir/OUTPUT_FILES/" --to_dir="$rundir/OUTPUT_FILES/"   \
       --from_template='${nt}.${sta}.BX${comp}.semd'  \
       --to_template='${nt}.${sta}.BX${comp}.sem.ascii'

echo "python $MEASURE_LIB/ascii2sac.py $current_dir/src_rec/STATIONS_${evtid}_globe  \
      $current_dir/src_rec/CMTSOLUTION_${evtid} $SYN_DIR"
\rm -f $SYN_DIR/*.sac # clean all sac file
python $MEASURE_LIB/ascii2sac.py $current_dir/src_rec/STATIONS_${evtid}_globe  \
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


# synthetic observations/STF
telemod_dir=$MEASURE_LIB/tele/
ncpu=8
echo " "
if [ $run_opt -ne 1 ]; then
  echo "estimating STF and synthetic data ..."
  mpirun -np $ncpu python $MEASURE_LIB/preprocess_tele.py $iter $evtid $run_opt
else
  echo "synthetic data ..."
  mpirun -np $ncpu python $MEASURE_LIB/preprocess_tele.py $iter $evtid $run_opt

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
  mkdir -p OUTPUT_FILES
  \cp ../../DATA/MEASUREMENT.PAR  .
  $MEASURE_LIB/bin/measure_adj
  coef=1.
  avgmap=`head -1 average_amp.dat`
  dt=`tail -1 average_amp.dat`
  coef=`echo $dt $avgmap |awk '{print $1/$2/$2}'`
  awk -v a=$evtid -v b=$coef '{$1=a;$29=$29*b; print}' window_chi >  $current_dir/misfits/$mod.${evtid}${lsflag}_${band}_tele_window_chi

  # continue for line search
  if [ $run_opt -eq 2 ]; then
    \rm -rf OUTPUT_FILES
    continue
  fi

  # convert sac to hdf5
  python $FWATLIB/measure/sac2h5.py ../seismogram.obs.$band.h5 *.sac.obs
  python $FWATLIB/measure/sac2h5.py ../seismogram.syn.$band.h5 *.sac.syn

  # rotate adjoint source 
  cd OUTPUT_FILES
  for f in `ls *iker*.adj |cut -d'.' -f1,2 |sort -n |uniq`;
  do 
    newf=`echo $f |awk -F'.' '{new_var=$2"."$1; print new_var}'`
    baz=`saclst baz f ../$newf.BXZ.sac.obs | awk '{print $2 * 3.1415926535898/180.}'`
    avgmap=`head -1 ../average_amp.dat`
    python $MEASURE_LIB/tele/cal_adjoint_src.py $f.BXZ.iker02.adj ../stf_Z.sac $avgmap
    python $MEASURE_LIB/tele/cal_adjoint_src.py $f.BXR.iker02.adj ../stf_R.sac $avgmap

    $MEASURE_LIB/bin/rotate_adj_src $baz $f.BXZ.iker02.adj $f.BXT.iker02.adj $f.BXR.iker02.adj $newf.BXE.adj $newf.BXN.adj 
    cp $f.BXZ.iker02.adj $newf.BXZ.adj
  done
  cd ../..
done

# continue for finish line search
if [ $run_opt -eq 2 ]; then
  exit 0;
fi

cd $current_dir
bash $MEASURE_LIB/sum_adj_source.sh $iter $evtid tele