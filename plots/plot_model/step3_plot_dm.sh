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

get_plotting_range() {
  local filename=$1
  local info=`awk '{print $1,$2,$4}' $filename | gmt gmtinfo  -C`
  local lonmin=`echo $info| awk '{print $1}'`
  local lonmax=`echo $info| awk '{print $2}'`
  local latmin=`echo $info| awk '{print $3}'`
  local latmax=`echo $info| awk '{print $4}'`
  local hmin=`echo $info| awk '{print $5}'`
  local hmax=`echo $info| awk '{print $6}'`

  local dlat=`echo "$latmin $latmax" | awk '{print ($2 - $1)}'`
  local dlon=`echo "$lonmin $lonmax" | awk '{print ($2 - $1)}'`

  local return_value="2 $latmin $latmax"
  if (( $(echo "$dlat < $dlon" | bc -l) )); then 
    return_value="1 $lonmin $lonmax"
  fi

  echo "$return_value $hmin $hmax"
}

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
    M=$(awk -v a="$vmin" -v b="$vmax" 'BEGIN{
        if (a < 0) a = -a;
        if (b < 0) b = -b;
        print (a > b ? a : b)
    }')
    vmin=$(awk -v M="$M" 'BEGIN{print -M}')
    vmax=$M
    echo $filename $vmin $vmax $vmin $vmax
    # if [ "$param" == "G0" ];  then 
    #   vmin=0
    #   vmax=0.05
    # fi
    gmt makecpt -T$vmin/$vmax/50+n -Z -D -Cvik -I > out.cpt
    #gmt grd2cpt $filename -Z -D -Cpolar -I  > out.cpt

    # plot
    proj=-JX12c/6c

    # check if we plot lat/lon 
    info=`get_plotting_range profiles/topo.verti.$ip.txt`
    rowid=`echo $info | awk '{print $1}'`
    xmin1=`echo $info | awk '{print $2}'`
    xmax1=`echo $info | awk '{print $3}'`
    ymin1=`echo $info | awk '{print $4}'`
    ymax1=`echo $info | awk '{print $5}'`
    bounds1=-R$xmin1/$xmax1/$ymin1/$ymax1
    if [ "$rowid" == "1" ]; then
      awk '{print $1,$4}' profiles/topo.verti.$ip.txt > temp.txt 
      tag="Longitude"
    else 
      awk '{print $2,$4}' profiles/topo.verti.$ip.txt > temp.txt 
      tag="Latitude"
    fi
    gmt begin pics/$param.diff.iter$idx.$name.$ip jpg 

      gmt basemap $bounds1 -JX12c/2c  -Bxaf+l"$tag" -Byaf+l"Elevation,m" -BWbrN -Y10c -X10c
      gmt plot temp.txt -W1p,black 

      gmt basemap $bounds $proj  -Bxaf+l"Distance,km" -Byaf+l"Depth,km" -BWSet -Y-6c
      gmt grdimage $filename -Cout.cpt -E200
      gmt colorbar -G$vmin/$vmax -Cout.cpt -Bxaf+l"$param"
    gmt end 

    \rm temp.txt
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
    M=$(awk -v a="$vmin" -v b="$vmax" 'BEGIN{
        if (a < 0) a = -a;
        if (b < 0) b = -b;
        print (a > b ? a : b)
    }')
    vmin=$(awk -v M="$M" 'BEGIN{print -M}')
    vmax=$M

    echo $filename $vmin $vmax $vmin $vmax
    gmt makecpt -T$vmin/$vmax/50+n -Z -D -Cvik -I > out.cpt
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