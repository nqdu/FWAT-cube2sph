#!/bin/bash
set -e

SHELL_HEADER_SEM(){
  if [[ "$PLATFORM"  == "local" ]];  then
    cat << EOF
#!/bin/bash
EOF
  elif [[ "$PLATFORM"  == "slurm" ]]; then 
  local narray=$1
  local stype=$2
  local walltime=$3
    cat << EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --array=1-$narray%10
#SBATCH --time=$walltime
#SBATCH --job-name=ADJ.$flag.$stype
#SBATCH --output=LOG/ADJ.$flag.$stype-%j_set%a.txt
#SBATCH --account=rrg-liuqy
#SBATCH --partition=compute
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
EOF
  elif [[ "$PLATFORM"  == "pbs" ]]; then 
    cat << EOF
#!/bin/bash -l
#PBS -l nodes=2:ppn=80
#PBS -l walltime=01:30:00
#PBS -N ADJ.$flag.$stype
#PBS -J 1-$narray%10
#PBS -q starq
#PBS -j oe
EOF
  else 
    echo "not implemented!"
  fi
}

SHELL_HEADER_POST() {
  if [[ "$PLATFORM"  == "local" ]];  then 
    cat << EOF
#!/bin/bash
EOF
  else
    cat << EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=00:15:59
#SBATCH --job-name=POST
#SBATCH --output=LOG/POST_%j.txt
#SBATCH --account=rrg-liuqy
#SBATCH --partition=compute
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
EOF
  fi 
}

SHELL_HEADER_WOLFE(){
  if [ "$PLATFORM"  == "local" ];  then 
    cat << EOF
#!/bin/bash
EOF
  else
    cat << EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=00:15:59
#SBATCH --job-name WOLFE
#SBATCH --output=LOG/WOLFE_%j.txt
#SBATCH --account=rrg-liuqy
#SBATCH --partition=compute
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
EOF
  fi 
}

RUN_SEM()
{
  # local vars
  local iter=$1
  local job_ids=()

  local SIMU_TYPES=(`fwat-utils getparam simulation/types|tr -d '[]",'\'`)
  
  # check size of njobs per job array == number of simu types
  if [[ "$PLATFORM"  != "local" ]];  then 
    if [ "${#NJOBS_PER_JOBARRAY[@]}" -ne "${#SIMU_TYPES[@]}" ]; then 
      echo "Error: The number of entries in NJOBS_PER_JOBARRAY (${#NJOBS_PER_JOBARRAY[@]}) does not match the number of simulation types (${#SIMU_TYPES[@]}). Please check your configuration."
      exit 1
    fi
  fi

  local nsimtypes="${#SIMU_TYPES[@]}"
  for ((isim=0;isim<$nsimtypes;isim++)); 
  do 
    local njobs=${NJOBS_PER_JOBARRAY[$isim]}
    local simu_type=${SIMU_TYPES[$isim]}
    local nevts=`awk 'END { print NR }' ${FWAT_SRC_REC}/sources.dat.$simu_type`
    local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`
    if [ "$PLATFORM"  == "local" ]; then 
      njobs=$nevts
    fi

    # copy files
    local fwd=tmp.adj.$simu_type.sh
    if [[ $simu_type == "noise" ]]; then
      SHELL_HEADER_SEM $narray $simu_type 00:25:00 > $fwd 
    else
      SHELL_HEADER_SEM $narray $simu_type 00:25:00 > $fwd 
    fi
    cat sbash_measure.sh >>  $fwd

    # substitute 
    \cp DATA/Par_file.$simu_type DATA/Par_file

    # run forward/adjoint simulation
    echo "forward/adjoint $simu_type simulation ..."

    # submit and get job id
    if [ "$PLATFORM"  == "local"  ]; then 
      bash $fwd $simu_type > LOG/ADJ.$flag.$simu_type.$iter.txt
    elif [ "$PLATFORM"  == "slurm" ]; then
      local jid=$(sbatch $fwd $simu_type |cut -d ' ' -f4 ) 
      job_ids+=($jid)
    else 
      local jid=$(qsub $fwd | cut -d '.' -f1)
      job_ids+=($jid)
    fi
  done 

  # return dependency strings if required
  if [[ "$PLATFORM"  != "local" ]];  then 
    local depend_string=$(IFS=:; echo "${job_ids[*]}")
    job_adj=$depend_string
  fi
}

RUN_POST ()
{
  local fwd=tmp.post.sh 
  SHELL_HEADER_POST 00:16:00 > $fwd 
  cat sbash_post.sh >> $fwd 

  echo "post processing ..."
  if [[ "$PLATFORM"  == "local" ]];  then  
    bash $fwd post > LOG/POST.$iter.txt 
  elif [[ "$PLATFORM"  == "slurm" ]]; then
    if [ "$flag" == "GRAD" ]; then 
      job_post=$(sbatch $fwd post | cut -d ' ' -f4)
    else 
      job_post=$(sbatch --dependency=afterok:${job_adj} $fwd  post | cut -d ' ' -f4)
    fi
  else 
    if [ "$flag" == "GRAD" ]; then 
      job_post=$(qsub $fwd | cut -d '.' -f1)
    else 
      job_post=$(qsub -W depend=afterok:$job_adj $fwd | cut -d '.' -f1)
    fi
  fi 
}

RUN_WOLFE () {

  local fwd=tmp.wolfe.sh 
  SHELL_HEADER_WOLFE 00:16:00 > $fwd 
  cat sbash_post.sh >> $fwd 

  # check wolfe condition
  echo "checking wolfe condition ..."

  if [[ "$PLATFORM"  == "local" ]];  then  
    bash $fwd wolfe > LOG/WOLFE.$iter.txt 
  elif [[ "$PLATFORM"  == "slurm" ]]; then
    job_post=$(sbatch --dependency=afterok:${job_adj} $fwd wolfe | cut -d ' ' -f4)
  else 
    job_post=$(qsub -W depend=afterok:$job_adj $fwd wolfe | cut -d '.' -f1)
  fi
} 

WAIT_FINISH() {
  if [[ "$PLATFORM"  == "local" ]];  then 
    echo ""
  elif [[ "$PLATFORM"  == "slurm" ]]; then 
    echo "waiting for post processing job $job_post to finish ..."
    while squeue -j $job_post | grep -q $job_post; do
      sleep 60
    done

    # the job may fail, check job status
    sleep 5
    # -X ignores job steps (like .batch), -n removes headers, -P removes whitespace formatting
    local final_state=$(sacct -j "$job_post" -X -n -P -o State | head -n 1)
    
    # Slurm states sometimes append text (e.g., "CANCELLED by 1234"). This strips it to just the first word.
    final_state=${final_state%% *}

    if [[ "$final_state" != "COMPLETED" ]]; then
      echo "Error: Slurm job $job_post failed with status: $final_state"
      exit 1
    fi
  else 
    echo "waiting for post processing job $job_post to finish ..."
    while qstat -u $USER | grep -q $job_post; do
      sleep 60
    done

    if ! qstat -x -f "$job_post" 2>/dev/null | grep -q "exit_status = 0"; then
      echo "Error: PBS job $job_post failed, was cancelled, or returned a non-zero exit status."
      exit 1
    fi

  fi
}
######### USER PARAMETERS ###############

source config.env
source utils.sh

# sanity check
SANITY_CHECK

# mkdir 
mkdir -p misfits optimize solver LOG

# some jobid 
job_adj=0
job_post=0

for ii in `seq 1 4`;do 

  # current model
  iter=`fwat-utils getparam iter ${LBFGS_FILE}`
  flag=`fwat-utils getparam flag ${LBFGS_FILE}`
  mod=M`printf %02d $iter`
  #mkdir -p LOG/${mod}
  echo "iteration $iter $mod $flag"

  # copy first model to MODEL_M00
  if [[ "$iter" -eq 0 && ! -d "${FWAT_OPT_DIR}/MODEL_M00" ]]; then
    mkdir -p ${FWAT_OPT_DIR}/MODEL_M00
    \cp initial_model/* ${FWAT_OPT_DIR}/MODEL_M00
  fi

  # check flag type and run 
  if [ $flag == "INIT" ]; then 
    RUN_SEM $iter
    
    # sum kernels, get search direction, generate trial model 
    RUN_POST

  elif [ $flag == "GRAD"  ];then 
    # get search direction, generate trial model 
    RUN_POST

  else  # line search
    RUN_SEM $iter

    RUN_WOLFE
  fi

  WAIT_FINISH
done