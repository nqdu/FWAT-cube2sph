#/bin/bash
set -e 
bounds=-R-122.0/-117.0/34.75/37.5 
proj=-JM12c 

# generate
lat0=35.9
lon0=-121.75
lat1=36.68
lon1=-117.15


gmt project -C$lon0/$lat0 -E$lon1/$lat1 -G1 -Q  > profile.dat 
dmax=`tail -1 profile.dat | awk '{print $3}' `

# interpolate
mkdir -p grdfolder pics 
if [ ! -f depth.dat  ];then
	awk '{print $3}' $file | sort -n | uniq > depth.dat
fi
depths=(`cat depth.dat`)
zmin=${depths[0]}
nz=${#depths[@]}
zmax=${depths[$nz-1]}

# plot vertical profile
bounds=-R0/$dmax/$zmin/$zmax
awk '{print $3,$5,$4}' out.dat | gmt surface -I200+n/200+n $bounds -Gout.grd 

vmin=`gmt grdinfo out.grd | grep v_min | awk '{print $3}' `
vmax=`gmt grdinfo out.grd | grep v_min | awk '{print $5}' `
gmt makecpt -T$vmin/$vmax/50+n -Z -Cseis > out.cpt
gmt begin out pdf
gmt basemap $bounds -JX12c -Bxa50f25 -Bya50f25 -BWSen
gmt grdimage out.grd   -Cout.cpt
gmt colorbar -Cout.cpt -Bxa500f250
gmt end 

rm out*.dat out.cpt out.grd 
