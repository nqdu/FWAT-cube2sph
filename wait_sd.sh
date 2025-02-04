#!/bin/bash


MODEL=M50

. parameters.sh 

icur=$(echo $MODEL |awk -F'M' '{print $2}')
inext=$(printf "%02d" `echo $MODEL |awk -F'M' '{print $2+1}'`)

  
rm -rf ./optimize/MODEL_M$inext 
mv ./optimize/MODEL_${MODEL}_step01 ./optimize/MODEL_M$inext 


