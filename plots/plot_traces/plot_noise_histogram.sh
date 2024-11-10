#!/bin/bash
module load gmt-6.0.0
set -e

if [ $# != 1 ];then 
    echo "./this Model(M00) "
    exit 1
fi 

mydir=../../solver/
proj=-JX12c
sc=0.5c/-1
M0=M25
M1=$1
seisdir=seismograms_comp
mkdir -p $seisdir

# get bands
fwat_file=../../fwat_params/FWAT.PAR.noise
SHORT_P=(`cat $fwat_file |grep 'SHORT_P:' |awk -F: '{print $2}'`)
LONG_P=(`cat $fwat_file |grep 'LONG_P:' |awk -F: '{print $2}'`)
NUM_FILTER=`echo ${#SHORT_P[@]}`
period_used=""
for ((i=0;i<$NUM_FILTER;i++));
do
  strt=`printf "T%03d_T%03d " ${SHORT_P[$i]} ${LONG_P[$i]}`
  period_used="$period_used""$strt"
done

for band in $period_used;
do 
  #residule file
  echo $band
  file0=./residual0.txt
  file1=./residual1.txt
  :> $file0
  :> $file1 
  cat ../../src_rec/sources.dat.noise |while read line;
  do 
    name=`echo $line |awk '{print $1}'`
    file=../../misfits/$M0.${name}_${band}_noise_window_chi 
    #MT and XCTT
    cat $file |awk '{if($13!=0) print $13}' >> $file0
    cat $file |awk '{if($13==0&&$15!=0) $15}' >> $file0

    file=../../misfits/$M1.${name}_${band}_noise_window_chi 
    #file=1.window_chi
    cat $file |awk '{if($13!=0) print $13}' >> $file1
    cat $file |awk '{if($13==0&&$15!=0) $15}' >> $file1 
  done

  # plot histogram
  bounds=-R-40/40/-1/100
  gmt begin $seisdir/dt_${band} jpg
  gmt basemap $bounds $proj -Bxaf+l"dT" -Byaf -BWSen
  gmt histogram $file0 -T0.4 -F -S -W1p,black
  gmt histogram $file1 -T0.4 -F -S -W1p,red
  gmt end 

  exit 1
done  