#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:59:59
#SBATCH --job-name FORWARD
#SBATCH --output=FORWARD_%j.txt
#SBATCH --partition=debug
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
#SBATCH --mail-type=ALL

set -e 

#####################################
###### Parameters ##################

# event sets
setb=12
sete=22

# simu type
simu_type=tele

# HPC facilities
nnodes=4
ntasks_per_node=40

# FKSEM path
fksem='/home/l/liuqy/nqdu/specfem3d/'

###### Parameters END ################
####################################

# LOAD all module
source activate base
module load NiaEnv/2019b
module load intel openmpi

echo "========================="
echo "Synthetic begin `date`"
echo "========================="

# copy parameters
if [ $simu_type == 'noise' ];then
  echo "copying params for noise"
  cp DATA/Par_file.noise DATA/Par_file
  cp fwat_params/FWAT.PAR.noise fwat_params/FWAT.PAR
  cp fwat_params/MEASUREMENT.PAR.noise fwat_params/MEASUREMENT.PAR
  cp src_rec/sources_ls.dat.noise src_rec/sources_ls.dat
elif [ $simu_type == 'tele' ];then
  echo "copying params for tele"
  cp DATA/Par_file.tele DATA/Par_file
  cp fwat_params/FWAT.PAR.tele fwat_params/FWAT.PAR
  cp fwat_params/MEASUREMENT.PAR.tele fwat_params/MEASUREMENT.PAR
  cp src_rec/sources_ls.dat.tele src_rec/sources_ls.dat
fi 

# machinefile 
NPROC=`grep "^NPROC" DATA/Par_file | awk '{print $3}'`
srun hostname | sort -n | uniq | awk -F. '{print $1}' > slurm.host 
python utils/generate_hostfile.py $nnodes $ntasks_per_node $NPROC
nodes=`ls slurm.host.? | wc -l`

# generate mesh if required
LOCAL_PATH="OUTPUT_FILES/DATABASES_MPI/"
mkdir -p $LOCAL_PATH
rm -rf $LOCAL_PATH/*.bin
./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH DATA/Par_file
./utils/change_par_file.sh LOCAL_PATH $LOCAL_PATH  DATA/meshfem3D_files/Mesh_Par_file
cp checkerboard_model_40km/* ${LOCAL_PATH}/ -r 
echo 'mpirun --oversubscribe -np $NPROC $fksem/bin/xmeshfem3D'
mpirun -machinefile slurm.host --oversubscribe -np $NPROC $fksem/bin/xmeshfem3D
echo 'mpirun --oversubscribe -np $NPROC $fksem/bin/xgenerate_databases'
mpirun -machinefile slurm.host --oversubscribe -np $NPROC $fksem/bin/xgenerate_databases

# run forward simulation
python src/forward.py M00 $setb $sete $simu_type $NPROC $nodes 

# add noise
# python src/add_noise_to_synthetics.py data 0.01
