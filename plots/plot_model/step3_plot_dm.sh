#!/bin/bash
set -e
source parameters.sh
source module_env_gmt

mkdir -p pics

if [ "$INTP_KL" == "1" ]; then 
  param_set=$KERNEL_SET
else 
  param_set=$MODEL_SET
fi

lsflag=""
if [ "$INTP_LS" == "1" ]; then 
  lsflag=".ls"
fi 

# vertical 
for param in $param_set; do 
for iter in $run_indx; do 
  if [ "$iter" == "0" ]; then 
    continue 
  fi 

  ii=`printf %02d $iter`
  idx=${ii}${lsflag}
  name=verti
  nfiles=`ls input/ |grep $name.*.loc |wc -l`
  for ip in `seq 1 $nfiles`; do
    filename=grdfiles/$param.diff.iter$idx.$name.$ip.grd
    grdc=grdfiles/$param.iter$idx.$name.$ip.grd
    grd0=grdfiles/$param.iter00${lsflag}.$name.$ip.grd
    gmt grdmath $grdc $grd0  SUB $grd0  DIV 100 MUL = $filename

    xmin=`gmt grdinfo $filename |grep x_min |awk '{print $3}'`
    xmax=`gmt grdinfo $filename |grep x_min |awk '{print $5}'`
    ymin=`gmt grdinfo $filename |grep y_min |awk '{print $3}'`
    ymax=`gmt grdinfo $filename |grep y_min |awk '{print $5}'`
    bounds=-R$xmin/$xmax/$ymin/$ymax

    info=`gmt grdinfo $filename -C`
    vmin=`echo $info| awk '{print $6}'`
    vmax=`echo $info| awk '{print $7}'`
    echo $filename $vmin $vmax $vmin $vmax
    # if [ "$param" == "G0" ];  then 
    #   vmin=0
    #   vmax=0.05
    # fi
    gmt makecpt -T$vmin/$vmax/50+n -Z -D -Cpolar -I > out.cpt
    #gmt grd2cpt $filename -Z -D -Cpolar -I  > out.cpt

    # plot
    proj=-JX12c/6c

    gmt begin pics/$param.diff.iter$idx.$name.$ip jpg 
      gmt basemap $bounds $proj  -Bxaf+l"Distance,km" -Byaf+l"Depth,km" -BWSet
      gmt grdimage $filename -Cout.cpt -E200
      gmt colorbar -G$vmin/$vmax -Cout.cpt -Bxaf+l"$param"
    gmt end 

  done 
done
done 


# vertical 
for param in $param_set; do 
for iter in $run_indx; do 
  if [ "$iter" == "0" ]; then 
    continue 
  fi 
  ii=`printf %02d $iter`
  idx=${ii}${lsflag}
  name=horiz
  nfiles=`ls input/ |grep $name.*.loc |wc -l`
  for ip in `seq 1 $nfiles`; do
    filename=grdfiles/$param.diff.iter$idx.$name.$ip.grd
    grdc=grdfiles/$param.iter$idx.$name.$ip.grd
    grd0=grdfiles/$param.iter00${lsflag}.$name.$ip.grd
    gmt grdmath $grdc $grd0  SUB $grd0  DIV 100 MUL = $filename
    
    xmin=`gmt grdinfo $filename |grep x_min |awk '{print $3}'`
    xmax=`gmt grdinfo $filename |grep x_min |awk '{print $5}'`
    ymin=`gmt grdinfo $filename |grep y_min |awk '{print $3}'`
    ymax=`gmt grdinfo $filename |grep y_min |awk '{print $5}'`
    bounds=-R$xmin/$xmax/$ymin/$ymax

    info=`gmt grdinfo $filename -C`
    vmin=`echo $info| awk '{print $6}'`
    vmax=`echo $info| awk '{print $7}'`
    # echo $filename $vmin $vmax
    # if [ "$param" == "G0" ];  then 
    #   vmin=0
    #   vmax=0.05
    # fi
    echo $filename $vmin $vmax $vmin $vmax
    gmt makecpt -T$vmin/$vmax/50+n -Z -D -Cpolar -I > out.cpt
    #gmt grd2cpt $filename -Z -D -Cpolar -I  > out.cpt

    # plot
    proj=-JM12c

    gmt begin pics/$param.diff.iter$idx.$name.$ip jpg
      gmt basemap $bounds $proj  -Bxaf+l"Distance,km" -Byaf+l"Depth,km" -BWSet
      gmt grdimage $filename -Cout.cpt -E200
      gmt coast -A200 -W1p,black
      gmt colorbar -G$vmin/$vmax -Cout.cpt -Bxaf+l"$param"
    gmt end 

  done 
done
done 