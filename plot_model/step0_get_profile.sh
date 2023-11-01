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
module load gcc openmpi
mkdir -p grdfolder pics 

# parameters
NPROC=160

# profile
lon0=-125; lon1=-119.5
#lon0=-124.25; lon1=-122
lat0=44.42; lat1=44.42
#lat0=45.5; lat1=45.5
z0=-220
z1=0.
dist=5
nz=101

# gmt grdcut @earth_relief_03s -R-125/-119.5/44/45  -Gout.grd
# horizontal line
python src/generate_line.py $lon0 $lon1 $lat0 $lat1 300 profile.txt

for name in A; do 
  :>input/prof$name.txt
  for ((i=0;i<$nz;i++));
  do
      dep=`printf %g $(echo "scale=4; 1000*($z0 + ($z1 - $z0) / ($nz-1) * $i)"|bc)`
      awk -v a=$dep '{print $1,$2,a,$3}' profile.txt >> input/prof$name.txt
  done

  # get grd file and plot
  awk '{print $1,$2,$3}' input/prof$name.txt > temp.txt
  mpirun -np $NPROC --oversubscribe ./bin/generate_slice ../../optimize/MODEL_M00 temp.txt false
  \rm temp.txt 
  mv temp.txt.loc input/prof$name.loc

done 
