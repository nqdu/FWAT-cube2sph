#!/bin/bash
#SBATCH --nodes=4 
#SBATCH --ntasks-per-node=20
#SBATCH --time=00:15:00
#SBATCH --job-name MESH
#SBATCH --output=MESH%j.txt
#SBATCH --partition=TH_HPC3
sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./OUTPUT_FILES/DATABASES_MPI" DATA/Par_file
sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./OUTPUT_FILES/DATABASES_MPI" DATA/meshfem3D_files/Mesh_Par_file

fksem='/THL8/home/iggluyf/nqdu/specfem3d-joint/'
NPROC=`grep "^NPROC" DATA/Par_file | awk '{print $3}'`
echo "mpirun --oversubscribe -np $NPROC $fksem/bin/xmeshfem3D"
mpirun -np $NPROC  $fksem/bin/xmeshfem3D
echo "mpirun --oversubscribe -np $NPROC $fksem/bin/xgenerate_databases"
mpirun -np $NPROC $fksem/bin/xgenerate_databases
