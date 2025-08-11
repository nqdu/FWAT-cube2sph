#!/bin/bash
set -e

SHELL_HEADER_SEM(){
  if [[ "$PLATFORM"  == "local" ]];  then
    cat << EOF
#!/bin/bash
EOF
  else
  local narray=$1
  local stype=$2
  local walltime=$3
    cat << EOF
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --array=1-$narray%5
#SBATCH --time=$walltime
#SBATCH --job-name=FWD.$stype
#SBATCH --output=LOG/FWD.$stype-%j_set%a.txt
#SBATCH --account=rrg-liuqy
#SBATCH --partition=compute
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
EOF
  fi
}
RUN_SEM()
{
  nsimtypes="${#SIMU_TYPES[@]}"
  for ((isim=0;isim<$nsimtypes;isim++)); 
  do 
    local njobs=${NJOBS_PER_JOBARRAY[$isim]}
    local simu_type=${SIMU_TYPES[$isim]}
    local nevts=`awk 'END { print NR }' src_rec/sources.dat.$simu_type`
    local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`
    if [ "$PLATFORM"  == "local" ]; then 
      njobs=$nevts
    fi

    # copy files
    fwd=tmp.fwd.$simu_type.sh
    if [[ $simu_type == "noise" ]]; then
      SHELL_HEADER_SEM $narray $simu_type 00:25:00 > $fwd 
    else
      SHELL_HEADER_SEM $narray $simu_type 00:25:00 > $fwd 
    fi
    cat sbash_forward.sh >>  $fwd

    # substitute 
    \cp DATA/Par_file.$simu_type DATA/Par_file

    # run forward/adjoint simulation
    echo "forward $simu_type simulation  ..."

    # submit and get job id
    if [ "$PLATFORM"  == "local"  ]; then 
      bash $fwd $simu_type > LOG/FWD.$simu_type.0.txt
    else 
      sbatch $fwd $simu_type
    fi
  done 
}

######### USER PARAMETERS ###############
source parameters.sh

# mkdir 
mkdir -p misfits optimize solver LOG

RUN_SEM