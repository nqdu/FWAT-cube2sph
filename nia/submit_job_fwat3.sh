#!/bin/bash
mod=$1
simu_type=$2
env=slurm # slurm or pbs
partition=debug
. checkpoint.sh 

###
if [ $# -ne 2 ];then
    echo " Usage: ./submit_fwat1.bash M?? noise/tele "
    echo "narg: " $#
    exit
fi
###


if [ $env == 'slurm' ];then
   fwd=sbash_fwat3_linesearch.sh 
elif [ $env == 'pbs' ];then
   fwd=pbs_fwat3_linesearch.sh 
fi

if [[ $simu_type == "tele" ]];then 
  cp src_rec/sources_ls.dat.tele src_rec/sources_ls -r 
else
  cp src_rec/sources_ls.dat.noise src_rec/sources_ls -r 
fi  


for ipart in `grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`;do
    #for ipart in 0.020;do
  echo "======================================="
  echo Model:$mod Part:$ipart
  echo "========================================"
  if [ $env == 'slurm' ];then
    #######   for slurm ######
    sed -i "/#SBATCH --job-name/c\#SBATCH --job-name LS${ipart}" $fwd 
    if [ $partition == "debug" ]; then 
      sed -i "/#SBATCH --time/c\#SBATCH --time=00:59:59" $fwd
    else 
      sed -i "/#SBATCH --time/c\#SBATCH --time=01:30:59" $fwd
    fi
    sed -i "/#SBATCH --output/c\#SBATCH --output=LS${ipart}_\%j.txt" $fwd
    sed -i "/MODEL=/c\MODEL=${mod}_step${ipart}" $fwd
    sed -i "/SET=/c\SET=ls" $fwd
    sed -i "/SIMU_TYPE=/c\SIMU_TYPE=${simu_type}" $fwd
    sed -i "/#SBATCH --partition/c\#SBATCH --partition=$partition" $fwd
    sbatch $fwd
    sleep 5

    if [[ $partition == "debug"  ]]; then
      CheckDebug $partition
    fi
  elif [ $env == 'pbs' ];then
    #######   for PBS ######
    sed -i "/#PBS -N/c\#PBS -N ${mod}.LS${ipart}" $fwd 
    sed -i "/#PBS -l walltime/c\#PBS -l walltime=01:30:00" $fwd
    sed -i "/MODEL=/c\MODEL=${mod}_step${ipart}" $fwd
    sed -i "/SET=/c\SET=ls" $fwd
    sed -i "/SIMU_TYPE=/c\SIMU_TYPE=${simu_type}" $fwd
    qsub $fwd
    sleep 3
  fi
done
##################
