#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:59:59
#SBATCH --job-name FWI
#SBATCH --output=FWI_%j.txt

set -e 

#####################################
###### Parameters ##################

if [ $# != 2  ];then
	echo "Usage ./this simu_type startmod(M00)"
	exit 1
fi

# simu_type
simu_type=$1

# start model
startmod=$2

###### Parameters END ################
####################################

# LOAD all module
#source activate base
#module load NiaEnv/2018a

# go with tele data 
startidx=$(printf %d `echo $startmod | cut -d'M' -f2`)
jj=`printf %02d $startidx`
# compute misfits for each step size
cd plots 
bash plot_misfit/plt_line_search.multiband.tele.bash M$jj

# next index
jnext=$(printf %02d $(echo "1+$startidx"|bc))

# find the min misfit and copy to Model_next
minval=`awk '{print $2}' M${jj}.mis.avg | sort -n | head -1`
step=`grep $minval M${jj}.mis.avg -n | cut -d: -f2 | awk '{print $1}'`
cd ../optimize
cp MODEL_M${jj}_step${step} MODEL_M$jnext -rf 
cd ..
echo "misfit for iteration $ii : $minval"

# save log file
mkdir -p LOGS/MODEL_$jj
mv  *.txt LOGS/MODEL_$jj -f

echo " "
