#!/bin/bash

set -e
source parameters.sh 
source module_env_gmt

if [ "$INTP_KL" == "1" ]; then 
  param_set=$KERNEL_SET
else 
  param_set=$MODEL_SET
fi

lsflag=""
if [ "$INTP_LS" == "1" ]; then 
  lsflag=".ls"
fi 

mkdir -p grdfiles

echo "$param_set txt to grd ..."

# make mask grd
gmt xyz2grd profiles/mask.dat -Ggrdfiles/mask.grd -I256+n/256+n -R$LON0_H/$LON1_H/$LAT0_H/$LAT1_H

# for vertical slice, get topography
bounds=-R$LON0_H/$LON1_H/$LAT0_H/$LAT1_H
gmt grdcut @earth_relief_30s -Ggrdfiles/topo.grd $bounds -Vq
for ip in `seq 1 $NSLICE_VERTI`; do 
  let ii=$p-1
  lon0=${LON0_V[$ii]}
  lon1=${LON1_V[$ii]}
  lat0=${LAT0_V[$ii]}
  lat1=${LAT1_V[$ii]}
  python src/generate_gc.py $lon0 $lon1 $lat0 $lat1 300 profile.txt
  gmt grdtrack profile.txt -Ggrdfiles/topo.grd > profiles/topo.verti.$ip.txt
done 


for param in $param_set ;do 
for iter in $run_indx;do 
  ii=`printf %02d $iter`
  idx=${ii}${lsflag}
  for name in verti horiz; do 
    nfiles=`ls input/ |grep $name.*.loc |wc -l`

    for ip in `seq 1 $nfiles`; do
      
      info=`gmt gmtinfo -C profiles/$param.iter$idx.$name.$ip.txt`
      x0=`echo $info | awk '{print $1}'`
      x1=`echo $info | awk '{print $2}'`
      z0=`echo $info | awk '{print $3}'`
      z1=`echo $info | awk '{print $4}'`
      dx=`echo "$x0 $x1" | awk '{print (-$1 + $2) / 255.}'`
      dz=`echo "$z0 $z1" | awk '{print (-$1 + $2)/255.}'`
      bounds=-R$x0/$x1/$z0/$z1 

      if [ "$INTP_KL" == "1" ]; then  
        
        vmin=`echo $info |awk '{print $5}'`
        vmax=`echo $info |awk '{print $6}'`
        flag=`echo "$vmin $vmax" |awk '{print $2-$1 >=0}'`
        #flag=$(echo "(-1*$vmin) >= $vmax" |bc)
        if [ "$flag" == "1" ]; then 
          vmax=`echo $vmin | awk '{print -$1}'`
        fi
        awk -v a=$vmax '{print $1,$2,$3/a}' profiles/$param.iter$idx.$name.$ip.txt > tmp.3 
      else 
        cat profiles/$param.iter$idx.$name.$ip.txt > tmp.3 
      fi

      gmt surface tmp.3 -Ggrdfiles/$param.iter$idx.$name.$ip.grd -I$dx/$dz $bounds -Vq

      if [ "$name" == "horiz" ]; then 
        gmt grdmath grdfiles/$param.iter$idx.$name.$ip.grd grdfiles/mask.grd MUL =  tmp.grd 
        \mv tmp.grd grdfiles/$param.iter$idx.$name.$ip.grd
      fi 
      
      \rm tmp.3
      gmt grdinfo -C grdfiles/$param.iter$idx.$name.$ip.grd
    done 
  done 
done 
done 