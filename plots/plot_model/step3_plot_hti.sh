#!/bin/bash
set -e
source module_env_gmt
source parameters.sh

#for material in hess_kernel;

for material in $MODEL_SET; 
do
  for iter in $run_indx; # model index
  do
    ii=`printf %02d $iter`
    for i in `seq 1 $NSLICE_HORIZ`; # profile index
    do
      # get model difference
      filename=grdfiles/$material.iter$ii.horiz.$i.grd

      xmin=`gmt grdinfo $filename |grep x_min |awk '{print $3}'`
      xmax=`gmt grdinfo $filename |grep x_min |awk '{print $5}'`
      ymin=`gmt grdinfo $filename |grep y_min |awk '{print $3}'`
      ymax=`gmt grdinfo $filename |grep y_min |awk '{print $5}'`
      bounds=-R$xmin/$xmax/$ymin/$ymax
      echo $bounds

      # compute g0 sintheta and g0 costheta
      paste profiles/G0.iter$ii.horiz.$i.txt profiles/phi.iter$ii.horiz.$i.txt  |  \
        awk '{print $1,$2,0.25*cos($6),0.25*sin($6)}' > aniso.txt.1
      cat aniso.txt.1| sort -n |awk 'NR % 4 == 0'  > aniso.txt
 
      #gmt grdmath grdfiles/G0.iter$ii.prof$prof.grd grdfiles/theta.iter$ii.prof$prof.grd

      # cut 
      gmt grdcut $filename -Gtmp.grd $bounds 
      filename=tmp.grd

      # makecpt
      info=`gmt grdinfo $filename -C`
      vmin=`echo $info| awk '{print $6}'`
      vmax=`echo $info| awk '{print $7}'`
      echo $vmin $vmax
      #vmin=0
      # vmin=-100
      # vmax=100
      gmt makecpt -T$vmin/$vmax/50+n -Z -D -Cvik -I > out.cpt
      gmt grd2cpt $filename -Z -D -Cvik -I  > out.cpt

      # plot
      proj=-JM12c

      gmt begin pics/${material}_iter${ii}.horiz.$i.hti  jpg

        # gmt basemap $bounds -JX12c/1.5c -Bxaf+l"Longitude,Deg" -Byaf+l"Topo,m" \
        #     -BWNeb+t"${prof}-${prof}1" -Y12c

        gmt basemap $bounds $proj -Bxaf -Byaf -BWSen+t"Depth = ${depth[$i]}km"
        
        gmt grdimage $filename -Cout.cpt -E200
        #gmt coast -A20 -W1p,black
        gmt velo aniso.txt $bounds $proj  -Sn0.2 -W0.5p,black

        #gmt plot station.lst $bounds $proj -St0.5c -Gred

        #awk '{print $1,$2/1000}' input/slab.bnd |gmt plot $bounds $proj -W1p,black
        #awk '{print $5,$3}' slab2.txt |gmt plot -W1.0p,black
        #gmt grdcontour $filename -C100
        gmt colorbar -G$vmin/$vmax -Cout.cpt -Bxaf+l"$material,km/s"
      gmt end

      \rm tmp.grd
    done
  done
done
