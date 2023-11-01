#!/bin/bash

source activate pygmt 

# profile
lon0=-124.25; lon1=-122
lat0=44.4; lat1=44.4
lonmin=-125; 
latmin=44
lonmax=-119.5
latmax=45
dist=5

# get topo 
bounds=-R$lonmin/$lonmax/$latmin/$latmax 

gmt grdcut @earth_relief_03s $bounds -Gout.grd 
:>profile.txt 
for i in `seq 0 499`;
do
  x=`echo "scale=6; $lon0 + ($lon1 - $lon0) / 499 * $i" |bc`
  echo $x $lat0 >> profile.txt 
done 
#gmt project -C$lon0/$lat0 -E$lon1/$lat1 -G$dist -Q > profile.txt 
awk '{print $1,$2}' profile.txt |gmt grdtrack  $bounds -Gout.grd > temp.txt  
awk '{print $1,$3}' temp.txt >  topo.txt 

info=`gmt gmtinfo topo.txt -C`
xmin=`echo $info | awk '{print $1}'`
xmax=`echo $info | awk '{print $2}'`
hmin=`echo $info | awk '{print $3}'`
hmax=`echo $info | awk '{print $4}'`
hmax=`echo "$hmax*1.1" |bc`

bounds=-R$xmin/$xmax/$hmin/$hmax
proj=-JX12c/2c
gmt begin out jpg
gmt basemap $bounds $proj -Bxaf+l"Distance,km" -Bya400f200+l"Topography,m"
gmt plot topo.txt -W0.5p,black
gmt end 
\rm out.grd profile.txt  temp.txt 