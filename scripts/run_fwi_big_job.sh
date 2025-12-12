#!/bin/bash

SHELL_HEADER_PRE() {
  local nodes=$1
  local nproc_per_node=$2 
  if [[ "$PLATFORM"  == "local" ]];  then 
    cat << EOF
#!/bin/bash
EOF
  elif [[ "$PLATFORM"  == "slurm" ]];  then 
    cat << EOF
#!/bin/bash
#SBATCH --nodes=$nodes
#SBATCH --ntasks-per-node=$nproc_per_node
#SBATCH --time=00:15:59
#SBATCH --job-name=FWI
#SBATCH --output=LOG/FWI_%j.log 
#SBATCH --account=rrg-liuqy
#SBATCH --partition=compute
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
EOF
  elif [[ "$PLATFORM"  == "pbs" ]];  then 
    cat << EOF
#!/bin/bash -l
#PBS -l nodes=${nodes}:ppn=$nproc_per_node
#PBS -l walltime=00:15:59
#PBS -N FWI
#PBS -q starq
#PBS -j oe
EOF
  else 
      echo "not implemented!"
  fi
}

source module_env 
. parameters.sh

# check input args 
# submit and get job id
if [ "$PLATFORM"  == "local"  ]; then 
  if [ "$#" -ne 2 ]; then
    echo "Usage: $0 sbash_measure NPROCS_TOTAL max_iter"
    exit 1
  fi
  NODES=1
  NPROCS_TOTAL=$1
  max_iter=$2
else 
  if [ "$#" -ne 3 ]; then
    echo "Usage: $0 sbash_measure NODES NPROCS_PER_NODE max_iter"
    exit 1
  fi
  NODES=$1
  NPROCS_PER_NODE=$2
  max_iter=$3
fi

# generate a temp file to submit job 
fwd=tmp.fwi.sh 
SHELL_HEADER_PRE $NODES $NPROCS_PER_NODE > $fwd
cat sbash_or_local_fwi.sh >> $fwd
chmod +x $fwd

exit 1

# submit job
if [ "$PLATFORM"  == "local"  ]; then 
  ./$fwd
elif [ "$PLATFORM"  == "slurm"  ]; then 
  sbatch $fwd
elif [ "$PLATFORM"  == "pbs"  ]; then 
  qsub $fwd
else 
  echo "not implemented!"
  exit 1
fi