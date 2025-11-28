#!bin/bash

# binary/python file location
SEM_PATH=~/specfem3d-cube2sph  # solver location
MPIRUN=mpirun

PLATFORM="local" # local/slurm

# simulation types
SIMU_TYPES=("noise")
SIMU_TYPES_USER_WEIGHT=(1.)
NJOBS_PER_JOBARRAY=(1)

############## STOP HERE ###################

change_par() 
{
  # get input args
  local param=$1
  local value=$2
  local file=$3

  # locate parameter
  oldstr=`grep "^$param " $file`
  newstr="$param           =     $value"

  sed  "s?$oldstr?$newstr?g" $file  > $file.temporary
  mv $file.temporary $file
}
