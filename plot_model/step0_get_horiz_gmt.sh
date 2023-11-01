#!/bin/bash
#SBATCH --nodes=4
#SBATCH --partition=debug
#SBATCH --ntasks=160
#SBATCH --time=00:59:59
#SBATCH --output=step0.log
set -e 
module load intel openmpi
source activate pygmt 
mkdir -p grdfolder pics 

function plot_horiz
{
    local isplot=$3
    local depth=$2
    local file=$1
    # convert txt to grd
    awk '{print $1,$2,$4/1000}' $file|  gmt surface -Ggrdfolder/$file.grd -I128+n/128+n $bounds 

    if $isplot ; then 
        # makecpt
        info=`gmt grdinfo grdfolder/$file.grd -C`
        local vmin=`echo $info| awk '{print $6}'`
        local vmax=`echo $info| awk '{print $7}'`
        #vmin=-10
        #vmax=10

        #echo $vmin $vmax
        echo $vmin $vmax 
        gmt makecpt -T$vmin/$vmax/50+n -Z -D -Cseis > out.$file.cpt 
        gmt begin pics/$file jpg 
        gmt basemap $bounds $proj -Bxa1f0.5 -Bya1f0.5 -BWSen+t"depth=${depth}"
        gmt grdimage grdfolder/$file.grd -Cout.$file.cpt -E200
        #gmt grdcontour grdfolder/$file.grd -C0.1 -A0.1 
        awk '{print $2,$1}' stations.dat | gmt plot -St0.5c -Gblack
        gmt coast -A500 -W1p,black
        gmt colorbar -Cout.$file.cpt  -Bxaf
        gmt end 
        rm out.$file.cpt 
    fi
}

# parameters
NPROC=160
param=vs 

# horizontal slice 
lonmin=-125; 
latmin=42
lonmax=-119.5
latmax=46.75
# lonmin=-122.0; 
# latmin=34.75
# lonmax=-117.0
# latmax=37.5

# get grd file and plot
bounds=-R$lonmin/$lonmax/$latmin/$latmax 
proj=-JM12c
for i in `seq 0 4`;
do 
    jj=`printf %02d $i`
    for dep in 30; do 
        python src/generate_horiz.py $lonmin $lonmax $latmin $latmax -$dep 100 100 input/$param.iter$jj.$dep.cords
        mpirun -np $NPROC --oversubscribe ./bin/read_media ../../optimize/MODEL_M$jj/ \
            input/$param.iter$jj.$dep.cords $param $param.iter$jj.dep$dep
        # mpirun -np $NPROC --oversubscribe ./bin/read_media ../tomo_models/model.checkerboard \
        #      input/$param.iter$jj.$dep.cords $param $param.iter$jj.dep$dep
        plot_horiz $param.iter$jj.dep$dep $dep true
    done 
     
done 


