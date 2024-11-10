#!/bin/bash

if [[ $# -ne 2  && $# -ne 3 ]]; then 
  echo "Usage: ./clean.sh MODEL evtid (deepclean)"
  exit 1
fi

mod=$1
evtid=$2
deepclean=0
if [ $# -eq 3 ]; then 
  deepclean=$3
fi 

# clean
rundir=solver/$mod/$evtid/
LOCAL_PATH=`grep '^LOCAL_PATH' $rundir/DATA/Par_file | awk '{print $3}'`
\rm -rf $rundir/$LOCAL_PATH
\rm -rf $rundir/OUTPUT_FILES/*.sac

# delete STF
fwat_file=$rundir/DATA/FWAT.PAR
SHORT_P=(`cat $fwat_file |grep 'SHORT_P:' |awk -F: '{print $2}'`)
LONG_P=(`cat $fwat_file |grep 'LONG_P:' |awk -F: '{print $2}'`)
NUM_FILTER=`echo ${#SHORT_P[@]}`
for ((i=0;i<$NUM_FILTER;i++));
do
  band=`printf "T%03d_T%03d" ${SHORT_P[$i]} ${LONG_P[$i]}`
  \rm -rf  $rundir/OUTPUT_FILES/$band/stf*
done

if [ $deepclean -eq 1 ]; then 
  \rm -rf $rundir
fi 