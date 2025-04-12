#!/bin/bash
#SBATCH --nodes=4 
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:15:00
#SBATCH --job-name MESH
#SBATCH --output=MESH%j.txt
#SBATCH --partition=debug
##SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
##SBATCH --mail-type=ALL

module load intel openmpi
#cd $SLURM_SUBMIT_DIR

#sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./OUTPUT_FILES/DATABASES_MPI" DATA/Par_file
#sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./OUTPUT_FILES/DATABASES_MPI" DATA/meshfem3D_files/Mesh_Par_file
./utils/change_par_file.sh LOCAL_PATH ./OUTPUT_FILES/DATABASES_MPI DATA/Par_file
./utils/change_par_file.sh LOCAL_PATH ./OUTPUT_FILES/DATABASES_MPI DATA/meshfem3D_files/Mesh_Par_file

fksem='/home/l/liuqy/nqdu/specfem3d/'
NPROC=`grep "^NPROC" DATA/Par_file | awk '{print $3}'`
echo "mpirun --oversubscribe -np $NPROC $fksem/bin/xmeshfem3D"
mpirun --oversubscribe -np $NPROC $fksem/bin/xmeshfem3D
echo "mpirun --oversubscribe -np $NPROC $fksem/bin/xgenerate_databases"
mpirun --oversubscribe -np $NPROC $fksem/bin/xgenerate_databases
