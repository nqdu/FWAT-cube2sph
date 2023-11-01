#!/bin/bash 

set -e

#####################################
###### Parameters ##################

# simu_type
iter=$1
simu_type=$2
step=$3

# go with tele data 
#startidx=$(printf %d `echo $startmod | cut -d'M' -f2`)
jj=$(printf %02d `echo "$iter + 0" | bc`)

# compute misfits for each step size
cd plots
out=`bash plot_misfit/cal_misfit.sh M$jj $simu_type $step`

echo $out 