#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:12:59
#SBATCH --job-name step0
#SBATCH --output=step0.txt 
#SBATCH --mem=0
#SBATCH --partition=debug

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
#module load NiaEnv/2019b
#module load intel openmpi
module load gmt-6.0.0
mkdir -p grdfolder pics profiles

# parameters
NPROC=160
param=vs

for name in A;  do
  # get grd file and plot
  info=`gmt gmtinfo -C input/prof${name}.txt`
  dmax=`echo $info | awk '{print $8}'`
  z0=`echo $info | awk '{print $5/1000}'`
  z1=`echo $info | awk '{print $6/1000}'`
  bounds=-R0/$dmax/$z0/$z1 
  proj=-JX12c/6c

  for iter in `seq 0 6`;
  do 
      jj=`printf %02d $iter`
      echo $jj 
      ./bin/interp3D_serial  ../../optimize/MODEL_M$jj input/prof$name.loc $param $param.out $NPROC
      #./bin/interp3D_serial  ../../../cascadia/tomo_models/model.checkerboard input/prof$name.loc $param $param.out $NPROC

      awk '{print $4}' input/prof$name.txt > tmp.1 
      awk '{print $3,$4}' $param.out > tmp.2 
      paste tmp.1 tmp.2 > tmp.3 

      # convert txt to grd
      gmt surface tmp.3 -Ggrdfolder/$param.iter$jj.prof$name.grd -I128+n/128+n $bounds -Vq
      mv tmp.3 profiles/$param.iter$jj.prof$name
      \rm tmp.1 tmp.2 $param.out 
      echo " " 
  done

done