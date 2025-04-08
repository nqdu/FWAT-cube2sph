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

. parameters.sh

set -e

# clean
rundir=solver/$mod/$evtid/
LOCAL_PATH=`grep '^LOCAL_PATH' $rundir/DATA/Par_file | awk '{print $3}'`
\rm -rf $rundir/$LOCAL_PATH
\rm -rf $rundir/OUTPUT_FILES/*.sac
\rm -rf $rundir/DATA/meshfem*
#\rm -rf $rundir/SEM
#\rm -rf $rundir/GRADIENT

# delete STF
for band in $rundir/OUTPUT_FILES/T[0-9]*_T[0-9]*;
do
  \rm -rf  $band
done
\rm -rf  $rundir/OUTPUT_FILES/timestamp*

if [ $deepclean -eq 1 ]; then 
  \rm -rf $rundir
fi 
