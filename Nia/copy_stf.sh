#!/bin/bash

for i in `seq 23 28`;
do
    cp solver/M00.set$i/$i/OUTPUT_FILES/stf_pca001_R.sac src_rec/STF_${i}.sac
done
