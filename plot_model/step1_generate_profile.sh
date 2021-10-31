#/bin/bash
set -e 
bounds=-R-122.0/-117.0/34.75/37.5 
proj=-JM12c 

# profile 
lat0=35.9
lon0=-121.75
lat1=36.68
lon1=-117.15

# two inputfiles diff = (v2 - v1)/v1 * 100
file1=vp.iter00 
file2=vp.iter06

gmt project -C$lon0/$lat0 -E$lon1/$lat1 -G1 -Q  > profile.dat 
dmax=`tail -1 profile.dat | awk '{print $3}' `

function ProfileInterp
{
    local depth=$1
    local i=$2
    local isplot=$5
    local file1=$3
    local file2=$4

    # compute difference between two grdfiles
    local file=grdfolder/$file1.$file2.$i.diff.grd
    if [ ! -f $file ];then
        gmt grdmath grdfolder/$i.$file2.grd grdfolder/$i.$file1.grd SUB grdfolder/$i.$file1.grd DIV 100 MUL = $file
    fi

    # interpolate 
    gmt grdtrack -G$file profile.dat >tmp.$i.dat 
    local npts=`cat tmp.$i.dat|wc -l`
    for j in `seq $npts`; do echo ${depths[$i]}; done > tmp1.$i.dat 
    paste tmp.$i.dat tmp1.$i.dat > out.$i.dat 

    if $isplot ; then 
        # makecpt
        local vmin=`gmt grdinfo $file | grep z_min | awk '{print $3}' `
        local vmax=`gmt grdinfo $file | grep z_min | awk '{print $5}' `
        #echo $vmin $vmax
        gmt makecpt -T$vmin/$vmax/50+n -Z -Cseis -D > out.$i.cpt 
        gmt begin pics/$file1.$file2.$i.diff pdf 
        gmt basemap $bounds $proj -Bxa1f0.5 -Bya1f0.5 -BWSen+t"depth=${depths[$i]}"
        gmt coast -A1000 -W1p,black
        gmt grdimage $file -Cout.$i.cpt
        awk '{print $1,$2}' profile.dat | gmt plot -W0.5p,black 
        gmt colorbar -Cout.$i.cpt  -Bxa2f1
        gmt end 
        rm out.$i.cpt 
    fi

    rm tmp1.$i.dat tmp.$i.dat
    echo ${depths[$i]}
}

# getinfo and interpolate 
if [ ! -f depth.dat  ];then
	awk '{print $3}' $file | sort -n | uniq > depth.dat
fi
depths=(`cat depth.dat`)
nz=${#depths[@]}
zmin=${depths[0]}
zmax=${depths[$nz-1]}
if true; then 
for ((i=0;i<$nz;i++));
do
    ProfileInterp $depths $i $file1 $file2 true  &
done 
wait 

# plot vertical profile
mkdir -p profiles
:>profiles/$file1.$file2.diff
for i in ./out.*.dat; do 
    cat $i >> profiles/$file1.$file2.diff
done
rm out*.dat
fi 

bounds1=-R0/$dmax/$zmin/$zmax 
awk '{print $3,$5,$4}' profiles/$file1.$file2.diff | gmt surface -I200+n/200+n $bounds1 -Gout.grd -Vq

vmin=`gmt grdinfo out.grd | grep v_min | awk '{print $3}'`
vmax=`gmt grdinfo out.grd | grep v_min | awk '{print $5}'`
vmin=-5
vmax=5
gmt makecpt -T$vmin/$vmax/50+n -Z -D -Cseis > out.cpt 
iter=`echo $file2 | cut -d'r' -f2`

gmt begin out.$iter pdf,jpg
gmt basemap $bounds1 -JX12c -Bxa50f25 -Bya50f25 -BWSen+t"Iteration $iter"
gmt grdimage out.grd -Cout.cpt 
gmt colorbar -Cout.cpt -Bxa2f1+l"@~\144@~ lnv (%) " -F
gmt end 

rm  out.cpt out.grd 
