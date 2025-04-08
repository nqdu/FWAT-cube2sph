#!/bin/bash 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
set -e 
source activate pygmt 

# parameters
workdir=../../
MODs=0 # start model
MODf=$1 # end model
setnum=$2  # source set 

if [[ $# != 2 ]]; then 
    echo "Usage: ./this modelindx(0,1,2 ...) sourceset(1,2,3 ...)"
    exit
fi 

# step1: mkdir for this pair to store plots 
mkdir -p sources_set$setnum 

# step2: get period used
fwat_par=$workdir/fwat_params/FWAT.PAR
nt=`grep NUM_FILTER $fwat_par |cut -d':' -f2`
Tlow=($(grep SHORT_P $fwat_par |cut -d':' -f2,`echo "$nt+1" |bc `))
Thigh=($(grep LONG_P $fwat_par |cut -d':' -f2,`echo "$nt+1" |bc `))

# step3: read stationlist and plot 
modsidx=`printf %02d $MODs`
modfidx=`printf %02d $MODf`
cat $workdir/src_rec/sources_set$setnum.dat| while read sline;
do
    SOURCENAME=`echo $sline |awk '{print $1}'`
    echo "ploting $SOURCENAME ..."

    # compute dist range
    dists=(`saclst dist f $workdir/fwat_data/$SOURCENAME/*.sac | awk '{print $2}' | gmt gmtinfo -C`)
    time=(``)
    mindist=`echo "${dists[0]} * 0.95" | bc `
    maxdist=`echo "${dists[1]} * 1.05" | bc `

    # plot by gmt sac 
    bounds=-R0/200/$mindist/$maxdist
    proj=-JX12c/12c

    for ((it=0;it<$nt;it++));do 
        mkdir -p sources_set$setnum/
        # read observed data and filter in this periods range 
        freqmin=`echo "scale=6; 1.0/${Thigh[$it]}"|bc`
        freqmax=`echo "scale=6; 1.0/${Tlow[$it]}"|bc`

        band=`printf "T%03d_T%03d"  ${Tlow[$it]} ${Thigh[$it]}`
        #ls $workdir/solver/M${modsidx}.set${setnum}/${SOURCENAME}/OUTPUT_FILES/*[obs,syn].sac.$band
        sac << EOF
        r $workdir/solver/M${modsidx}.set${setnum}/${SOURCENAME}/OUTPUT_FILES/*[obs,syn].sac.$band
        ch LCALDA true
        w append .M${modsidx}.bp 
        q
EOF
        # move to dir
        mv $workdir/solver/M${modsidx}.set${setnum}/$SOURCENAME/OUTPUT_FILES/*.sac.$band.M${modsidx}.bp  sources_set$setnum 

        sac << EOF
        r $workdir/solver/M${modfidx}.set${setnum}/${SOURCENAME}/OUTPUT_FILES/*[obs,syn].sac.$band
        ch LCALDA true
        w append .M${modfidx}.bp 
        q
EOF
        mv $workdir/solver/M${modfidx}.set${setnum}/$SOURCENAME/OUTPUT_FILES/*.sac.$band.M${modfidx}.bp  sources_set$setnum 
        
        # add sac header 
        ./add_sacheader.sh $setnum $SOURCENAME $workdir

        # add 
        gmt begin sources_set$setnum/$SOURCENAME.${Tlow[$it]}_${Thigh[$it]} jpg 
        gmt basemap -BWSne+t"MODEL$modsidx" $bounds $proj -Bxaf+l"Time,s" -Byaf+l"Distance,km"
        gmt sac sources_set$setnum/*obs.sac.$band.M${modsidx}.bp -Ek -M1.5c -W0.5p,black $bounds $proj 
        gmt sac sources_set$setnum/*syn.sac.$band.M${modsidx}.bp -Ek -M1.5c -W0.5p,red $bounds $proj 

        gmt basemap -BWSne+t"MODEL$modfidx" $bounds $proj -Bxaf+l"Time,s" -Byaf+l"Distance,km" -X14c 
        gmt sac sources_set$setnum/*obs.sac.$band.M${modfidx}.bp -Ek -M1.5c -W0.5p,black $bounds $proj 
        gmt sac sources_set$setnum/*syn.sac.$band.M${modfidx}.bp -Ek -M1.5c -W0.5p,red $bounds $proj 
        gmt end 

        rm sources_set$setnum/*.bp
    done 
done 