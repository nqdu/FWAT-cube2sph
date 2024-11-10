#!/bin/bash
set_fwat1()
{
  local simu_type=$1
  local njobs=$2
  local start_set=$3
  local nevts=`cat src_rec/sources.dat.$simu_type|wc -l`
  local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`

  if [[ $iter == 0 || $iter == $iter_end ]]; then
    sed -i "/SHOW_DETAILS:/c\SHOW_DETAILS: .true." fwat_params/FWAT.PAR
  fi

  # substitute 
  cp sbash_fwat1_fwd_measure_adj.sh tmp.fwat1.$simu_type.sh
  \cp DATA/Par_file.$simu_type DATA/Par_file
  local fwd=tmp.fwat1.$simu_type.sh
  sed -i "/#SBATCH --array=/c\#SBATCH --array=1-$narray%5" $fwd
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/NJOBS=/c\NJOBS=$njobs" $fwd
  sed -i "/START_SET=/c\START_SET=$start_set" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd
  sed -i "/#SBATCH --output=/c#SBATCH --output=FWD_ADJ_${simu_type}-%j_set%a.txt" $fwd 

  if [[ $simu_type == "noise" ]]; then
    sed -i "/#SBATCH --time=/c\#SBATCH --time=01:00:00" $fwd
  else
    sed -i "/#SBATCH --time=/c\#SBATCH --time=00:25:00" $fwd
  fi

  # run forward/adjoint simulation
  echo "computing forward/adjoint simulation for $simu_type ..."
}

#!/bin/bash
set_fwat3()
{
  local simu_type=$1
  local njobs=$2
  local start_set=$3
  local nevts=`cat src_rec/sources.dat.ls.$simu_type|wc -l`

  # njobs *2
  if [ $njobs -lt 5 ]; then
    njobs=`echo "$njobs *2"|bc`
  fi
  local narray=`echo "($nevts + $njobs - 1) / $njobs"|bc`

  # substitute 
  cp sbash_fwat3_linesearch.sh tmp.fwat3.$simu_type.sh
  \cp DATA/Par_file.$simu_type DATA/Par_file
  local fwd=tmp.fwat3.$simu_type.sh
  sed -i "/#SBATCH --array=/c\#SBATCH --array=1-$narray%5" $fwd
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/NJOBS=/c\NJOBS=$njobs" $fwd
  sed -i "/START_SET=/c\START_SET=$start_set" $fwd
  sed -i "/simu_type=/c\simu_type=$simu_type" $fwd
  sed -i "/#SBATCH --output=/c#SBATCH --output=LS_${simu_type}-%j_set%a.txt" $fwd 

  if [[ $simu_type == "noise" ]]; then
    sed -i "/#SBATCH --time=/c\#SBATCH --time=01:00:00" $fwd
  else
    sed -i "/#SBATCH --time=/c\#SBATCH --time=00:25:00" $fwd
  fi

  # run forward/adjoint simulation
  echo "computing line search for $simu_type ..."
}

# parameters
iter_start=27
iter_end=27
FIRST_ITER=25

# L-BFGS params
lbfgs_start=25

# simu_type
simu_type1=tele
simu_type2=noise
NJOBS1=1
NJOBS2=4

# mkdir 
mkdir -p misfits optimize solver

# some jobid 
job_adj=0
job_post=0
job_step=0
job_line=0
job_wait=0

# set number
# initialize set number
setb=1

for iter in `seq $iter_start $iter_end`;do
  # current model
  mod=M`printf %02d $iter`
  mod_lbfgs_start=M`printf %02d $lbfgs_start`
  echo "iteration $iter $mod $mod_lbfgs_start"

  # create misfit file 
  :> plots/$mod.mis

  # run RF/tele simulation
  sete=`cat src_rec/sources.dat.$simu_type1 |wc -l`
  sete1=`echo "$setb + $sete - 1" |bc`
  set_fwat1 $simu_type1 $NJOBS1 $setb
  #job_adj1=$(sbatch tmp.fwat1.$simu_type1.sh|cut -d ' ' -f4 )

  sete=`cat src_rec/sources.dat.$simu_type2 |wc -l`
  sete2=`echo "$setb + $sete - 1" |bc`
  set_fwat1 $simu_type2 $NJOBS2 $setb
  #job_adj2=$(sbatch tmp.fwat1.$simu_type2.sh|cut -d ' ' -f4 )

  # run line search
  if [ ${mod} == $mod_lbfgs_start  ]; then 
    echo ${mod} > lbfgs.in
    echo "-1" >> lbfgs.in
  else 
    info=`head -1 lbfgs.in`
    echo $info > lbfgs.in 
    echo "1" >> lbfgs.in
  fi
  #exit

  # run postprocessing
  # merge source files
  :> ./src_rec/sources.dat.joint
  cat ./src_rec/sources.dat.$simu_type1 > ./src_rec/sources.dat.joint
  cat ./src_rec/sources.dat.$simu_type2 >> ./src_rec/sources.dat.joint
  \cp sbash_fwat2_postproc_opt.sh   tmp.fwat2.sh 
  fwd=tmp.fwat2.sh 
  sed -i "/MODEL=/c\MODEL=${mod}" $fwd
  sed -i "/SOURCE_FILE=/c\SOURCE_FILE=./src_rec/sources.dat.joint" $fwd

  # save misfit for iter 0
  if [[ "$iter" == "$FIRST_ITER" ]];then
    cat > tmp.dat << EOF
:> misfit0.log
:> optimize/weight_kl.txt
mis1=\$(python cal_misfit.py $mod $simu_type1 |cut -d' ' -f1)
mis2=\$(python cal_misfit.py $mod $simu_type2 |cut -d' ' -f1)
for iset in \$(seq 1 $sete1); do echo "1.0"|bc -l >> optimize/weight_kl.txt; done
for iset in \$(seq 1 $sete2); do echo "\$mis1 / \$mis2"|bc -l >> optimize/weight_kl.txt; done 
EOF
    # substitute a line in $fwd
    nline=`grep -n ^PRECOND $fwd |cut -d: -f1`
    let nline=nline-1
    sed -n "1,${nline}p" $fwd > tmp.fwat2.cal.sh 
    cat tmp.dat  >> tmp.fwat2.cal.sh
    \rm tmp.dat
    let nline=nline+1
    sed -n "${nline},\$p" $fwd >> tmp.fwat2.cal.sh 
    \mv tmp.fwat2.cal.sh $fwd 
    #\rm tmp.fwat2.cal.sh
  fi
  echo "postprocessing ..."
  #job_post=$(sbatch --dependency=afterok:${job_adj1},${job_adj2} $fwd | cut -d ' ' -f4)

  # run line search
  set_fwat3 $simu_type1 $NJOBS1 $setb
  set_fwat3 $simu_type2 $NJOBS2 $setb
  #job_line=$(sbatch --dependency=afterok:${job_post} tmp.fwat3.$simu_type.sh | cut -d ' ' -f4)

  # generate next model
  echo "generating opt model ..."
  \cp sbash_fwat4_opt_model.sh tmp.fwat4.sh
  fwd=tmp.fwat4.sh
  modfirst=M`printf %02d $FIRST_ITER`
  cat > tmp.dat << EOF
mis1_0=\$(python cal_misfit.py $modfirst $simu_type1 00 |cut -d' ' -f1)
mis2_0=\$(python cal_misfit.py $modfirst $simu_type2 00 |cut -d' ' -f1)
mis1=\$(python cal_misfit.py $mod $simu_type1 00 |cut -d' ' -f1)
mis2=\$(python cal_misfit.py $mod $simu_type2 00 |cut -d' ' -f1)
chi=\$(echo "\$mis1/\$mis1_0 + \$mis2 / \$mis2_0" |bc -l)
mis1=\$(python cal_misfit.py $mod $simu_type1 01 |cut -d' ' -f1)
mis2=\$(python cal_misfit.py $mod $simu_type2 01 |cut -d' ' -f1)
chi1=\$(echo "\$mis1/\$mis1_0 + \$mis2 / \$mis2_0" |bc -l)
chi=\$(printf %f \$chi)
chi1=\$(printf %f \$chi1)
EOF
  # substitute a line in $fwd
  nline=`grep -n ^info= $fwd |cut -d: -f1|head -1`
  nline1=`grep -n ^chi1= $fwd |cut -d: -f1|head -1`
  let nline=nline-1
  sed -n "1,${nline}p" $fwd > tmp.fwat4.cal.sh 
  cat tmp.dat  >> tmp.fwat4.cal.sh
  \rm tmp.dat
  let nline1=nline1+1
  sed -n "${nline1},\$p" $fwd >> tmp.fwat4.cal.sh 
  \mv tmp.fwat4.cal.sh $fwd 
  #\rm tmp.fwat2.cal.sh
  sed -i "/model=/c\model=${mod}" $fwd
  #job_step=$(sbatch --dependency=afterok:${job_line} $fwd | cut -d ' ' -f4)
  exit 

  # wait job to finish
  srun --dependency=afterok:${job_step} --partition=compute --nodes=1 --ntasks=1 --time=00:15:02 wait.sh 
  mkdir -p LOG/$mod
  mv *.txt LOG/$mod
done
