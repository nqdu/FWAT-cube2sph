#!/bin/bash
set -e
source activate pygmt

for material in vs;
do
  for iter in `seq 11 11`; # model index
  do
    ii=`printf %02d $iter`
    echo $ii
    for prof in A; # profile index
    do
      # get model difference
      filename=grdfolder/$material.iter$ii.prof$prof.grd

      # makecpt
      info=`gmt grdinfo $filename -C`
      vmin=`echo $info| awk '{print $6}'`
      vmax=`echo $info| awk '{print $7}'`
      vmax=4800
      gmt makecpt -T$vmin/$vmax/50+n -Z -D -Cseis > out.cpt

      # plot
      xmin=`gmt grdinfo $filename |grep x_min |awk '{print $3}'`
      xmax=`gmt grdinfo $filename |grep x_min |awk '{print $5}'`
      ymin=`gmt grdinfo $filename |grep y_min |awk '{print $3}'`
      ymax=`gmt grdinfo $filename |grep y_min |awk '{print $5}'`
      bounds=-R$xmin/$xmax/$ymin/$ymax
      proj=-JX12c/6c

      # plot topo
      info=`gmt gmtinfo topo.txt -C`
      xmin=`echo $info | awk '{print $1}'`
      xmax=`echo $info | awk '{print $2}'`
      hmin=`echo $info | awk '{print $3}'`
      hmax=`echo $info | awk '{print $4}'`
      hmax=`echo "$hmax*1.1" |bc`
      bounds1=-R$xmin/$xmax/$hmin/$hmax


      gmt begin model_${material}_${ii}  jpg

        gmt basemap $bounds1 -JX12c/1.5c -Bxaf+l"Longitude,Deg" -Byaf+l"Topo,m" \
            -BWNeb+t"${prof}-${prof}1" -Y12c
        gmt plot topo.txt -W0.5p,black

        gmt basemap $bounds $proj -Bxaf+l"Distance,km" -Byaf+l"Depth,km" -Y-6.2c -BWSet
        gmt grdimage $filename -Cout.cpt -E200
        awk '{print $5,$3}' slab2.txt |gmt plot -W1.0p,black
        #gmt grdcontour $filename -C100
        gmt colorbar -G$vmin/$vmax -Cout.cpt -Bxaf+l"$material,km/s"
      gmt end
    done
  done
done