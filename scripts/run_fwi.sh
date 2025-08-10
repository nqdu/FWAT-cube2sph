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
#SBATCH --job-name=FWD_ADJ.$flag.$stype
#SBATCH --output=LOG/FWD_ADJ.$flag.$stype-%j_set%a.txt
#SBATCH --account=rrg-liuqy
#SBATCH --partition=compute
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
EOF
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
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:15:59
#SBATCH --job-name=POST
#SBATCH --output=LOG/POST_%j.txt
#SBATCH --account=rrg-liuqy
#SBATCH --partition=compute
#SBATCH --mem=12G
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
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:15:59
#SBATCH --job-name WOLFE
#SBATCH --output=LOG/WOLFE_%j.txt
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
  # local vars
  local iter=$1
  local job_ids=()

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
      bash $fwd $simu_type > LOG/FWD_ADJ.$flag.$simu_type.$iter.txt
    else 
      local jid=$(sbatch $fwd $simu_type |cut -d ' ' -f4 ) 
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
  cat sbash_postproc_kl.sh >> $fwd 

  echo "post processing ..."
  if [[ "$PLATFORM"  == "local" ]];  then  
    bash $fwd > LOG/POST.$iter.txt 
  else 
    job_post=$(sbatch --dependency=afterok:${job_adj} $fwd | cut -d ' ' -f4)
  fi 
}

RUN_WOLFE () {

  local fwd=tmp.wolfe.sh 
  SHELL_HEADER_WOLFE 00:16:00 > $fwd 
  cat sbash_wolfe.sh >> $fwd 

  # check wolfe condition
  echo "checking wolfe condition ..."

  if [[ "$PLATFORM"  == "local" ]];  then  
    bash $fwd > LOG/WOLFE.$iter.txt 
  else 
    job_post=$(sbatch --dependency=afterok:${job_adj} $fwd | cut -d ' ' -f4)
  fi
} 

WAIT_FINISH() {
  if [[ "$PLATFORM"  == "local" ]];  then 
    echo ""
  else 
    srun --dependency=afterok:${job_post} --nodes=1 --time=00:01:05 --ntasks=1 --job-name=wait  ./wait.sh
  fi
}

######### USER PARAMETERS ###############
source parameters.sh

# mkdir 
mkdir -p misfits optimize solver LOG

# some jobid 
job_adj=0
job_post=0

for ii in `seq 1 4`;do 

  # current model
  iter=`fwat-utils getparam iter fwat_params/lbfgs.yaml`
  flag=`fwat-utils getparam flag fwat_params/lbfgs.yaml`
  mod=M`printf %02d $iter`
  #mkdir -p LOG/${mod}
  echo "iteration $iter $mod $flag"

  # copy first model to MODEL_M00
  if [[ "$iter" -eq 0 && ! -d "optimize/MODEL_M00" ]]; then
    mkdir -p optimize/MODEL_M00
    \cp initial_model/* optimize/MODEL_M00
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