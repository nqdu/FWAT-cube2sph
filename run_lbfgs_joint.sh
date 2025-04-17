#!/bin/bash
set_fwat1()
{
  local simu_type=$1
  local njobs=$2
  local start_set=$3
  local nevts=`cat src_rec/sources.dat.$simu_type|wc -l`
  local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`

  # if [[ $iter == 0 || $iter == $iter_end ]]; then
  #   sed -i "/SHOW_DETAILS:/c\SHOW_DETAILS: .true." fwat_params/FWAT.PAR
  # fi

  # substitute 
  cp sbash_measure.sh tmp.fwat1.$simu_type.sh
  \cp DATA/Par_file.$simu_type DATA/Par_file
  local fwd=tmp.fwat1.$simu_type.sh
  sed -i "/#SBATCH --array=/c\#SBATCH --array=1-$narray%5" $fwd
  sed -i "/START_SET=/c\START_SET=$start_set" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd

  if [[ $simu_type == "noise" ]]; then
    sed -i "/#SBATCH --time=/c\#SBATCH --time=01:00:00" $fwd
  else
    sed -i "/#SBATCH --time=/c\#SBATCH --time=00:25:00" $fwd
  fi

  # run forward/adjoint simulation
  echo "forward/adjoint $simu_type simulation for new model  ..."
}

set_fwat2() 
{
	:> ./src_rec/sources.dat.joint
	cat ./src_rec/sources.dat.$simu_type1 > ./src_rec/sources.dat.joint
	cat ./src_rec/sources.dat.$simu_type2 >> ./src_rec/sources.dat.joint
	local fwd=tmp.fwat2.sh
  \cp sbash_postproc_kl.sh  $fwd 
  sed -i "/SOURCE_FILE=/c\SOURCE_FILE=./src_rec/sources.dat.joint" $fwd

  # save misfit for iter 0
  if [[ "$iter" == "$FIRST_ITER" ]];then
    local sete=`cat src_rec/sources.dat.$simu_type1 |wc -l`
    local sete1=`echo "$setb + $sete - 1" |bc`
    local sete=`cat src_rec/sources.dat.$simu_type2 |wc -l`
    local sete2=`echo "$setb + $sete - 1" |bc`
    cat > tmp.dat << EOF
:> misfit0.log
:> optimize/weight_kl.txt
mis1=\$(python \$MEASURE_LIB/cal_misfit.py $mod $simu_type1 |cut -d' ' -f1)
mis2=\$(python \$MEASURE_LIB/cal_misfit.py $mod $simu_type2 |cut -d' ' -f1)
for iset in \$(seq 1 $sete1); do echo "1.0"|bc -l >> optimize/weight_kl.txt; done
for iset in \$(seq 1 $sete2); do echo "\$mis1 / \$mis2"|bc -l >> optimize/weight_kl.txt; done 
EOF
    # substitute a line in $fwd
    local nline=`grep -n ^PRECOND $fwd |cut -d: -f1`
    let nline=nline-1
    sed -n "1,${nline}p" $fwd > tmp.fwat2.cal.sh 
    echo " " >> tmp.fwat2.cal.sh
    cat tmp.dat  >> tmp.fwat2.cal.sh
    echo " " >> tmp.fwat2.cal.sh
    \rm tmp.dat
    let nline=nline+1
    sed -n "${nline},\$p" $fwd >> tmp.fwat2.cal.sh 
    \mv tmp.fwat2.cal.sh $fwd 
    #\rm tmp.fwat2.cal.sh
  fi

  echo "post processing for new model  ..."
}

#!/bin/bash
set_fwat3()
{
  local simu_type=$1
  local njobs=$2
  local start_set=$3
  local nevts=`cat src_rec/sources.dat.$simu_type|wc -l`
  local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`

  # substitute 
  cp sbash_measure.sh tmp.fwat3.$simu_type.sh
  \cp DATA/Par_file.$simu_type DATA/Par_file
  local fwd=tmp.fwat3.$simu_type.sh
  sed -i "/#SBATCH --array=/c\#SBATCH --array=1-$narray%5" $fwd
  sed -i "/#SBATCH --job-name=/c\#SBATCH --job-name=LS" $fwd
  sed -i "/#SBATCH --output=/c\#SBATCH --output=LS-%j_set%a.txt" $fwd
  sed -i "/NJOBS=/c\NJOBS=$njobs" $fwd
  sed -i "/START_SET=/c\START_SET=$start_set" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd


  if [[ $simu_type == "noise" ]]; then
    sed -i "/#SBATCH --time=/c\#SBATCH --time=01:00:00" $fwd
  else
    sed -i "/#SBATCH --time=/c\#SBATCH --time=00:25:00" $fwd
  fi

  # run forward/adjoint simulation
  echo "forward/adjoint $simu_type simulation for new model  ..."
}

set_fwat4 () {
   echo "generating opt model ..."
  \cp sbash_wolfe.sh tmp.fwat4.sh
  local fwd=tmp.fwat4.sh
  local modfirst=M`printf %02d $FIRST_ITER`
  cat > tmp.dat << EOF
mis1_0=\$(python \$MEASURE_LIB/cal_misfit.py $modfirst $simu_type1 00 |cut -d' ' -f1)
mis2_0=\$(python \$MEASURE_LIB/cal_misfit.py $modfirst $simu_type2 00 |cut -d' ' -f1)
mis1=\$(python \$MEASURE_LIB/cal_misfit.py $mod $simu_type1 00 |cut -d' ' -f1)
mis2=\$(python \$MEASURE_LIB/cal_misfit.py $mod $simu_type2 00 |cut -d' ' -f1)
chi=\$(echo "\$mis1/\$mis1_0 + \$mis2 / \$mis2_0" |bc -l)
mis1=\$(python \$MEASURE_LIB/cal_misfit.py $mod $simu_type1 01 |cut -d' ' -f1)
mis2=\$(python \$MEASURE_LIB/cal_misfit.py $mod $simu_type2 01 |cut -d' ' -f1)
chi1=\$(echo "\$mis1/\$mis1_0 + \$mis2 / \$mis2_0" |bc -l)
chi=\$(printf %f \$chi)
chi1=\$(printf %f \$chi1)
EOF

  # substitute a line in $fwd
  local nline=`grep -n ^info= $fwd |cut -d: -f1|head -1`
  local nline1=`grep -n ^chi1= $fwd |cut -d: -f1|head -1`
  let nline=nline-1
  sed -n "1,${nline}p" $fwd > tmp.fwat4.cal.sh 
  cat tmp.dat  >> tmp.fwat4.cal.sh
  \rm tmp.dat
  let nline1=nline1+1
  sed -n "${nline1},\$p" $fwd >> tmp.fwat4.cal.sh 
  \mv tmp.fwat4.cal.sh $fwd 
  sed -i "/SIMU_TYPE=/c\SIMU_TYPE=joint" $fwd
}

set -e 
. parameters.sh

# simu_type
# simu_type
simu_type1=tele
simu_type2=noise
NJOBS1=1
NJOBS2=2

# initialize set number
setb=1
FIRST_ITER=0

########### STOP HERE ###################

# mkdir 
mkdir -p misfits optimize solver

# some jobid 
job_adj1=0
job_adj2=0
job_post=0
job_step=0
job_line=0
job_wait=0

# set number
for ii in `seq 1 1`;do 

  # current model
  iter=`python $FWATLIB/get_param.py iter $FWATPARAM/lbfgs.yaml`
  flag=`python $FWATLIB/get_param.py flag $FWATPARAM/lbfgs.yaml`
  mod=M`printf %02d $iter`
  echo "iteration $iter $mod $flag"

  # copy first model to MODEL_M00
  if [[ "$iter" -eq 0 && ! -d "optimize/MODEL_M00" ]]; then
    # set model name
    mkdir -p optimize/MODEL_M00
    \cp ./DATABASES_MPI/* optimize/MODEL_M00
  fi

  # check flag type and run 
  if [ $flag == "INIT" ]; then 
    set_fwat1 $simu_type1 $NJOBS1 $setb
    set_fwat1 $simu_type2 $NJOBS2 $setb
    job_adj1=$(sbatch tmp.fwat1.$simu_type1.sh|cut -d ' ' -f4 )
    job_adj2=$(sbatch tmp.fwat1.$simu_type2.sh|cut -d ' ' -f4 )
    
    # sum kernels, get search direction, generate trial model 
    set_fwat2
    fwd=tmp.fwat2.sh
    job_post=$(sbatch --dependency=afterok:${job_adj1},${job_adj2} $fwd | cut -d ' ' -f4)
    
  elif [ $flag == "GRAD"  ];then 
    # get search direction, generate trial model 
    set_fwat2
    fwd=tmp.fwat2.sh
    job_post=$(sbatch $fwd | cut -d ' ' -f4)

  else  # line search
    set_fwat3 $simu_type1 $NJOBS1 $setb
    set_fwat3 $simu_type2 $NJOBS2 $setb 
    job_adj1=$(sbatch tmp.fwat3.$simu_type1.sh|cut -d ' ' -f4 )
    job_adj2=$(sbatch tmp.fwat3.$simu_type2.sh|cut -d ' ' -f4 )

    # check wolfe condition
    set_fwat4
    fwd=tmp.fwat4.sh
    job_post=$(sbatch --dependency=afterok:${job_adj1},${job_adj2} $fwd | cut -d ' ' -f4)
  fi

  # wait to finish
  srun --dependency=afterok:${job_post} --time=00:00:10 --nodes=1 --time=00:00:10 --ntasks=1 --job-name=wait --partition=compute ./wait.sh
done
