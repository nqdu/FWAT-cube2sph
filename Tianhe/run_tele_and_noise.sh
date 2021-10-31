#!/bin/bash

# run noise
./submit_job_fwat1.sh M00 noise 1 10

# check if tasks are finished
idx=`squeue -u iggluyf | wc -l`
set -e
while [ $idx -ne 1  ]
do
	let j=$idx-1
	echo "$j tasks remaining ..."
	sleep 30
	idx=`squeue -u iggluyf | wc -l`
done

./submit_job_fwat1.sh M00 tele 12 22 
