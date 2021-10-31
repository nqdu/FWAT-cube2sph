#!/bin/bash
set -e 
bounds=-R-122.0/-117.0/34.75/37.5 
proj=-JM12c 

function gmtsurface 
{
    local depth=$1
    local i=$2
    local isplot=$4
    local file=$3

    # convert txt to grd 
    awk -v a=${depths[$i]} '$3==a {print $1,$2,$4}' $file |
    gmt surface -Ggrdfolder/$i.$file.grd -I100+n/100+n $bounds -Vq

    if $isplot ; then 
        # makecpt
        local vmin=`gmt grdinfo grdfolder/$i.$file.grd | grep v_min | awk '{print $3}' `
        local vmax=`gmt grdinfo grdfolder/$i.$file.grd| grep v_min | awk '{print $5}' `
        vmin=-5
        vmax=5

        #echo $vmin $vmax
        gmt makecpt -T$vmin/$vmax/50+n -Z -Cseis > out.$i.cpt 
        gmt begin pics/$$i.$file pdf 
        gmt basemap $bounds $proj -Bxa1f0.5 -Bya1f0.5 -BWSen+t"depth=${depths[$i]}"
        gmt coast -A1000 -W1p,black
        gmt grdimage grdfolder/$i.$file.grd -Cout.$i.cpt
        gmt colorbar -Cout.$i.cpt  -Bxa500f250
        gmt end 
        rm out.$i.cpt 
    fi

    echo ${depths[$i]}
}

# interpolate
mkdir -p grdfolder pics 
if [ ! -f depth.dat  ];then
	awk '{print $3}' vp.iter00 | sort -n | uniq > depth.dat
fi
depths=(`cat depth.dat`)
nz=${#depths[@]}

for file in vs.iter* ;do 
    echo $file
    for ((i=0;i<$nz;i++));
    do
        #echo $file
		gmtsurface $depths $i $file false & 
    done 
    wait 
done 
