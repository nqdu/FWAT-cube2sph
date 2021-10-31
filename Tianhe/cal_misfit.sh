#!/bin/bash

if [ $# -ne 2 ];then 
    echo "Usage ./this mod noise/tele"
    exit 1
fi 

modnum=$1
SIMU_TYPE=$2

# compute misfits for each step size
j=`printf %02d $modnum`
cd plots 
if [ $SIMU_TYPE == "noise" ];then 
    bash plot_misfit/plt_line_search.multiband.ANAT.bash M$j 
else 
    bash plot_misfit/plt_line_search.multiband.tele.bash M$j 
fi 

# next index
jnext=$(printf %02d $(echo "$modnum+1"|bc))

# find the min misfit and copy to Model_next
minval=`awk '{print $2}' M${j}.mis.avg | sort -n | head -1`
step=`grep $minval M${j}.mis.avg -n | cut -d: -f2 | awk '{print $1}'`
cd ../optimize
cp MODEL_M${j}_step${step} MODEL_M$jnext -rf 
cd ..
echo "misfit for iteration $ii : $minval"

# save log file
mkdir -p LOGS/MODEL_$j
mv output_fwat* FWD* POST*  LS* LOGS/MODEL_$j -f

echo " "