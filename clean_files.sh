#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:02:59
#SBATCH --job-name=CLEAN_FILES
#SBATCH --output=CLEAN_%j.txt
#SBATCH --partition=compute 

set -e

if [[ $# != 1  ]];then
    echo "Usage: ./clean_files num"
    exit 1
fi

i=$1

#for i in `seq 1 34`;
#do
  jj=`printf "%02d" $i`
  for f in solver/M$jj/*/*GRAD*;
  do 
    echo  $f 
   rm -rf $f & 
  done

  for f in solver/M$jj/*/SEM
  do
     echo  $f 
     rm -rf $f &
  done

  for f in solver/M${jj}.ls/*;
  do
    echo $f
    rm -rf $f &
  done 
#done 
