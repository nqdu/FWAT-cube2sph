#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:12:59
#SBATCH --job-name step0
#SBATCH --output=step0.txt 
#SBATCH --mem=0
#SBATCH --partition=debug

source parameters.sh
source module_env

rm -rf input
mkdir -p input 

# horizontal 
echo "" 
echo "working on horizontal profiles ..."
python src/generate_plane.py $LON0_H $LON1_H $LAT0_H $LAT1_H 256 profile.txt  
for i in `seq 1 $NSLICE_HORIZ`;
do 
  :> input/horiz.$i.txt  
  echo "horiz slice number $i"
  # 
  let ii=$i-1
  dep=${DEPTH_H[$ii]}
  awk  '{print "XZ ADF",$2,$1,0.,a*1000}' a=$dep profile.txt > temp.txt
  $cube2sph_dir/bin/write_stations_file temp.txt temp1.txt rotation_nu .false. .false.

  # merge data
  awk  '{print $4,$3,$6/1000}' temp.txt >  horiz.$i.txt.temp
  awk '{print $4,$3,$6}' temp1.txt >  temp.txt 
  paste temp.txt horiz.$i.txt.temp > input/horiz.$i.txt  

  # interpolate
  mpirun -np $NPROC  $specfem_dir/bin/xcreate_slice_loc input/horiz.$i.txt  $DATABASE_DIR temp.txt.loc
  \rm temp1.txt rotation_nu
  mv temp.txt.loc input/horiz.$i.loc
  \rm -f temp* horiz.$i.txt.temp 

done
\rm profile.txt


# vertical
echo "" 
echo "working on vertical profiles ..."
for i in `seq 1 $NSLICE_VERTI`;
do 

  echo "vertical slice number $i"

  let ii=$i-1
  lon0=${LON0_V[$ii]}
  lon1=${LON1_V[$ii]}
  lat0=${LAT0_V[$ii]}
  lat1=${LAT1_V[$ii]}
  python src/generate_gc.py $lon0 $lon1 $lat0 $lat1 300 profile.txt

  :>input/verti.$i.txt 
  :> verti.$i.temp
  nz=200
  for ((j=0;j<$nz;j++));
  do
    dep=`echo "$MAX_DEP $nz $j" |awk '{printf "%g", 1000*($1+(0-$1)/($2-1)*$3)}'`
    #dep=`printf %g $(echo "scale=6; 1000*($MAX_DEP + (0. - $MAX_DEP ) / ($nz-1.) * $j)"|bc)`
    awk  '{print "XZ ADF",$1,$2,a,$3}' a=$dep profile.txt >> verti.$i.temp
  done

  # get grd file and plot
  awk '{print $1,$2,$4,$3,0.0,-$5}' verti.$i.temp > temp.txt
  $cube2sph_dir/bin/write_stations_file temp.txt temp1.txt rotation_nu .false. .false.
  awk '{print $4,$3,$6}' temp1.txt > temp2.txt
  awk  '{print $5,$6}' verti.$i.temp  > temp3.txt 
  paste temp2.txt temp3.txt > input/verti.$i.txt 

  mpirun -np $NPROC  $specfem_dir/bin/xcreate_slice_loc input/verti.$i.txt  $DATABASE_DIR  temp.txt.loc

  \rm temp1.txt rotation_nu
  mv temp.txt.loc input/verti.$i.loc
  \rm temp* verti.$i.temp
  \rm profile.txt
done
