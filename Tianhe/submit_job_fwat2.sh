#!/bin/bash
mod=$1
setb=$2
sete=$3
env=slurm # slurm or pbs
partition=TH_HPC3

###
if [ $# -ne 3 ];then
    echo " Usage: ./submit_job_fwat2.bash M?? setb sete "
    echo "narg: " $#
    exit
fi
echo $mod
###

if [ $env == 'slurm' ];then 
  script=sbash_fwat2_postproc_opt.sh
elif [ $env == 'pbs' ];then
  script=pbs_fwat2_postproc_opt.sh
fi

if [ $mod == "M00" ];then
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./OUTPUT_FILES/DATABASES_MPI" DATA/Par_file
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./OUTPUT_FILES/DATABASES_MPI" DATA/meshfem3D_files/Mesh_Par_file
else
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${mod}" DATA/Par_file
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${mod}" DATA/meshfem3D_files/Mesh_Par_file
fi 


sed -i "/MODEL=/c\MODEL=${mod}" $script
sed -i "/SETB=/c\SETB=set${setb}" $script
sed -i "/SETE=/c\SETE=set${sete}" $script

if [ $env == 'slurm' ];then 
    sed -i "/#SBATCH --job-name/c\#SBATCH --job-name=POST" $script
    sed -i "/#SBATCH --time/c\#SBATCH --time=01:59:59" $script
    sed -i "/#SBATCH --output/c\#SBATCH --output=POST_\%j.txt" $script
    sed -i "/#SBATCH --partition/c\#SBATCH --partition=$partition" $script
    sbatch $script
elif [ $env == 'pbs' ];then
    qsub $script
fi
##################
