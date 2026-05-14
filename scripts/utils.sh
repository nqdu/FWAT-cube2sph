#!bin/bash

# functions used below 
####################################

# change parameter in file
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

# sanity check for inversion
SANITY_CHECK() 
{

  if [ "$PLATFORM" != "local" ] && [ "$PLATFORM" != "slurm" ] && [ "$PLATFORM" != "pbs" ]; then 
    echo "Error: Unsupported platform specified: $PLATFORM. Please set PLATFORM to 'local', 'slurm', or 'pbs' in parameters.sh."
    exit 1
  fi

  # check if SEM_PATH exists
  if [ ! -d "$SEM_PATH" ]; then
    echo "Error: SEM_PATH directory does not exist: $SEM_PATH. Please check your configuration."
    exit 1
  fi

  # check if fwat-utils is available
  if ! command -v fwat-utils &> /dev/null; then
    echo "Error: fwat-utils command not found. Please ensure it is installed and in your PATH."
    exit 1 
  fi

  # sanity check for simu types and weights
  local SIMU_TYPES=(`fwat-utils getparam simulation/types|  tr -d '[]",'\'`)
  local SIMU_WEIGHTS=(`fwat-utils getparam simulation/weights|  tr -d '[]",'\'`)
  
  # if more than 1 simulation type, check if norm_type is valid
  if [ "${#SIMU_TYPES[@]}" -gt 1 ]; then 
    local norm_type=`fwat-utils getparam simulation/norm_type`
    if [ "$norm_type" != "gradient" ] && [ "$norm_type" != "misfit" ]; then 
      echo "Error: Invalid norm_type specified: $norm_type. Please set norm_type to 'gradient' or 'misfit' when using multiple simulation types."
      exit 1
    fi
  fi

  if [ "${#SIMU_TYPES[@]}" -ne "${#SIMU_WEIGHTS[@]}" ]; then 
    echo "Error: The number of simulation types (${#SIMU_TYPES[@]}) does not match the number of weights (${#SIMU_WEIGHTS[@]}). Please check your configuration."
    exit 1
  fi
}