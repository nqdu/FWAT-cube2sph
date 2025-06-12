#!/bin/bash

set -e
source parameters.sh 
source module_env_gmt


lsflag=""
if [ "$INTP_LS" == "1" ]; then 
  lsflag=".ls"
fi 

mkdir -p grdfiles

param_set=$KERNEL_SET
echo "$param_set txt to grd ..."

for param in $param_set ;do 
for iter in $run_indx;do 
for name in verti horiz; do 
  ii=`printf %02d $iter`
  idx=${ii}${lsflag}
  nfiles=`ls input/ |grep $name.*.loc |wc -l`
  for evt in `ls $SOLVER_DIR/M$ii/`; do 
  for ip in `seq 1 $nfiles`; do
    info=`gmt gmtinfo -C profiles/$param.$evt.iter$idx.$name.$ip.txt`
    x0=`echo $info | awk '{print $1}'`
    x1=`echo $info | awk '{print $2}'`
    z0=`echo $info | awk '{print $3}'`
    z1=`echo $info | awk '{print $4}'`
    dx=`echo "$x0 $x1" | awk '{print (-$1 + $2) / 255.}'`
    dz=`echo "$z0 $z1" | awk '{print (-$1 + $2)/255.}'`
    bounds=-R$x0/$x1/$z0/$z1 
    vmin=`echo $info |awk '{print $5}'`
    vmax=`echo $info |awk '{print $6}'`
    flag=`echo "$vmin $vmax" |awk '{print $2-$1 >=0}'`
    #flag=$(echo "(-1*$vmin) >= $vmax" |bc)
    if [ "$flag" == "1" ]; then 
    vmax=`echo $vmin | awk '{print -$1}'`
    fi
    awk -v a=$vmax '{print $1,$2,$3/a}' profiles/$param.$evt.iter$idx.$name.$ip.txt > tmp.3 

    gmt surface tmp.3 -Ggrdfiles/$param.$evt.iter$idx.$name.$ip.grd -I$dx/$dz $bounds 
    \rm tmp.3
    gmt grdinfo -C grdfiles/$param.$evt.iter$idx.$name.$ip.grd
  done 
  done 

  # get info
  

done 
done 
done 