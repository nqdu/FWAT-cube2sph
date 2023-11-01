#!/bin/bash
set -e 

function change_par(){
  local param=$1
  local value=$2

  # locate parameter
  local oldstr=`grep "^$param" run_fwi.sh`
  local newstr="$param=$value"

  sed -i "s?$oldstr?$newstr?g" run_fwi.sh 
}

# parameters for set
setbt=23 # set tele start
setet=28  # set tele end
setbn=1 # set noise start
seten=22  # set noise end

# start and end model
iter_start=21
iter_end=30

# flag
iflag_tele=0

for ii in `seq $iter_start $iter_end`;
do
  echo "iteration" $ii $iflag_tele 
  if [[ $iflag_tele == 0 ]]; then 
    simu_type=noise
    setb=$setbn 
    sete=$seten 
  else 
    simu_type=tele 
    setb=$setbt 
    sete=$setet 
  fi

  # substitute variables 
  change_par simu_type $simu_type
  change_par setb $setb 
  change_par sete $sete 
  change_par iter_start $ii
  change_par iter_end $ii

  # run fwi
  bash run_fwi.sh 

  # change iflag_tele
  iflag_tele=`echo "($iflag_tele + 1) %2" |bc `

done 