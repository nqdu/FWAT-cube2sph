#!/bin/bash
# =====================================================================
#  run.sh  (monitor-integrated)
#
#  Same structure as the original run.sh, with:
#    * SHELL_HEADER_* functions kept verbatim
#    * RUN_SEM / RUN_POST / RUN_WOLFE rewritten to register jobs with
#      monitor.py so failures can be detected and resubmitted
#    * WAIT_FINISH delegated to monitor.py's wait loop
#
# =====================================================================
set -e

# ---------------------------------------------------------------------
#  Original header builders — kept as-is
# ---------------------------------------------------------------------
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

# ---------------------------------------------------------------------
#  Helper: look up the current scheduler ID of a logical job name.
#  (monitor.py rewrites this ID on resubmission, so always read fresh.)
# ---------------------------------------------------------------------
_jobid_of() {
  fwat-utils monitor --cmd=show | awk -v n="$1" '$1==n {print $2}'
}

# ---------------------------------------------------------------------
#  RUN_SEM — submits one job array per simulation type
# ---------------------------------------------------------------------
RUN_SEM()
{
  local iter=$1
  local job_names=()

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

    # build script
    local fwd=tmp.adj.$simu_type.sh
    if [  "$simu_type" == "noise" ]; then
      SHELL_HEADER_SEM $narray $simu_type 00:45:00 > $fwd
    else
      SHELL_HEADER_SEM $narray $simu_type 00:25:00 > $fwd
    fi
    cat sbash_measure.sh >>  $fwd

    # substitute
    \cp DATA/Par_file.$simu_type DATA/Par_file

    echo "forward/adjoint $simu_type simulation ..."

    # submit
    if [ "$PLATFORM"  == "local"  ]; then
      bash $fwd $simu_type > LOG/ADJ.$flag.$simu_type.$iter.txt
      continue
    fi

    local jid
    if [ "$PLATFORM"  == "slurm" ]; then
      jid=$(sbatch --parsable $fwd $simu_type)
    else
      jid=$(qsub -F "$simu_type" $fwd | cut -d '.' -f1 | tr -d '[]')
    fi

    # register with monitor so it can detect failures / resubmit
    local name="sem_${simu_type}_iter${iter}"
    fwat-utils monitor --cmd=register \
        --name        "$name" \
        --jobid       "$jid" \
        --script      "$fwd" \
        --array-spec  "1-$narray" \
        --script-args "$simu_type"

    job_names+=("$name")
  done

  # comma-separated list of logical parent names for downstream functions
  if [[ "$PLATFORM"  != "local" ]];  then
    job_adj_names=$(IFS=,; echo "${job_names[*]}")
  fi
}

# ---------------------------------------------------------------------
#  RUN_POST — depends on all SEM jobs from this iteration
# ---------------------------------------------------------------------
RUN_POST ()
{
  local fwd=tmp.post.sh
  SHELL_HEADER_POST 00:16:00 > $fwd
  cat sbash_post.sh >> $fwd

  echo "post processing ..."

  if [[ "$PLATFORM"  == "local" ]];  then
    bash $fwd post > LOG/POST.$iter.txt
    return
  fi

  # Resolve current scheduler IDs of parent SEM jobs (skip for GRAD)
  local -a parents_arg=()
  local -a dep_ids=()
  if [ "$flag" != "GRAD" ] && [ -n "$job_adj_names" ]; then
    IFS=',' read -ra parents_arg <<< "$job_adj_names"
    for n in "${parents_arg[@]}"; do
      dep_ids+=("$(_jobid_of "$n")")
    done
  fi

  local jid
  if [[ "$PLATFORM" == "slurm" ]]; then
    if [ ${#dep_ids[@]} -gt 0 ]; then
      local depstr=$(IFS=:; echo "${dep_ids[*]}")
      jid=$(sbatch --parsable --dependency=afterok:$depstr $fwd post)
    else
      jid=$(sbatch --parsable $fwd post)
    fi
  else  # pbs
    if [ ${#dep_ids[@]} -gt 0 ]; then
      local depstr=$(IFS=:; echo "${dep_ids[*]}")
      jid=$(qsub -W depend=afterok:$depstr -F post $fwd | cut -d '.' -f1)
    else
      jid=$(qsub -F post $fwd | cut -d '.' -f1)
    fi
  fi

  local post_name="post_iter${iter}"
  fwat-utils monitor --cmd=register \
      --name        "$post_name" \
      --jobid       "$jid" \
      --script      "$fwd" \
      --script-args "post" \
      ${parents_arg[@]:+--parents "${parents_arg[@]}"}

  job_post_name="$post_name"
}

# ---------------------------------------------------------------------
#  RUN_WOLFE — line-search check, depends on SEM jobs
# ---------------------------------------------------------------------
RUN_WOLFE () {
  local fwd=tmp.wolfe.sh
  SHELL_HEADER_WOLFE 00:16:00 > $fwd
  cat sbash_post.sh >> $fwd

  echo "checking wolfe condition ..."

  if [[ "$PLATFORM"  == "local" ]];  then
    bash $fwd wolfe > LOG/WOLFE.$iter.txt
    return
  fi

  local -a parents_arg=()
  local -a dep_ids=()
  IFS=',' read -ra parents_arg <<< "$job_adj_names"
  for n in "${parents_arg[@]}"; do
    dep_ids+=("$(_jobid_of "$n")")
  done
  local depstr=$(IFS=:; echo "${dep_ids[*]}")

  local jid
  if [[ "$PLATFORM" == "slurm" ]]; then
    jid=$(sbatch --parsable --dependency=afterok:$depstr $fwd wolfe)
  else
    jid=$(qsub -W depend=afterok:$depstr -F wolfe $fwd | cut -d '.' -f1)
  fi

  local wolfe_name="wolfe_iter${iter}"
  fwat-utils monitor --cmd=register \
      --name        "$wolfe_name" \
      --jobid       "$jid" \
      --script      "$fwd" \
      --script-args "wolfe" \
      --parents     "${parents_arg[@]}"

  job_post_name="$wolfe_name"
}

# ---------------------------------------------------------------------
#  WAIT_FINISH — hand off to the monitor's polling loop
# ---------------------------------------------------------------------
WAIT_FINISH() {
  if [[ "$PLATFORM"  == "local" ]];  then
    return 0
  fi

  echo "monitor: waiting for all registered jobs to finish ..."
  fwat-utils monitor --cmd=wait
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "monitor reported failure (exit $rc); aborting pipeline"
    exit $rc
  fi

  # Clear state for the next iteration so retry counters reset.
  fwat-utils monitor --cmd=reset
}

######### USER PARAMETERS ###############

source config.env
source utils.sh

# sanity check
SANITY_CHECK

# mkdir
mkdir -p misfits optimize solver LOG

# Reset monitor state at the start of a fresh pipeline run.
[[ "$PLATFORM" != "local" ]] && fwat-utils monitor --cmd=reset

# logical-name tracking variables used between RUN_* calls
job_adj_names=""
job_post_name=""

for ii in `seq 1 4`;do

  # current model
  iter=`fwat-utils getparam iter ${LBFGS_FILE}`
  flag=`fwat-utils getparam flag ${LBFGS_FILE}`
  mod=M`printf %02d $iter`
  echo "iteration $iter $mod $flag"

  # copy first model to MODEL_M00
  if [[ "$iter" -eq 0 && ! -d "${FWAT_OPT_DIR}/MODEL_M00" ]]; then
    mkdir -p ${FWAT_OPT_DIR}/MODEL_M00
    \cp initial_model/* ${FWAT_OPT_DIR}/MODEL_M00
  fi

  # check flag type and run
  if [ $flag == "INIT" ]; then
    RUN_SEM $iter
    echo " "

    # sum kernels, get search direction, generate trial model
    RUN_POST

  elif [ $flag == "GRAD"  ];then
    # get search direction, generate trial model
    RUN_POST

  else  # line search
    RUN_SEM $iter
    echo " "

    RUN_WOLFE
  fi

  echo ""

  WAIT_FINISH
done