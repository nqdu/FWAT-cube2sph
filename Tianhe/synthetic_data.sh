#/bin/bash
set -e 
. checkpoint.sh 

# clean files
LOCAL_PATH="OUTPUT_FILES/DATABASES_MPI/"
rm -rf $LOCAL_PATH/*.bin

# mesh and database
NPROC=`grep "^NPROC" DATA/Par_file | awk '{print $3}'`
echo "Generating database on $NPROC cores"
cp checkerboard_model_40km/* ${LOCAL_PATH}/ -r 
sbatch sbash_fwat0_mesh.sh

# check finished
fini=`ls | grep MESH |wc -l`
while [ $fini -ne 1 ]
do 
	sleep 5 
	fini=`ls | grep MESH |wc -l`
done 
sleep 30
rm MESH*

echo "Forward computing ..."
# synthetic data for checkerboard
#./submit_job_fwat_forward.sh M00 noise 1 10
#CheckForward M00 1 10

./submit_job_fwat_forward.sh M00 tele 12 22
CheckForward M00 12 22

python add_noise_to_synthetics.py data 0.0
