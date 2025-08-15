#!/bin/bash
specfem_dir="${HOME}/specfem3d-cube2sph/"
cube2sph_dir=$specfem_dir/utils/cube2sph/
NPROC=80
MPIRUN=mpirun

# horizontal slice
NSLICE_HORIZ=1
DEPTH_H=(225)  # add more if you want 
LON0_H=114; LON1_H=132;
LAT0_H=41; LAT1_H=46

# vertical slice
NSLICE_VERTI=1
MAX_DEP=400

# start points
LON0_V=(115 )
LAT0_V=(44)

#end points
LON1_V=(132 )
LAT1_V=(44)

# database_dir
DATABASE_DIR=../../DATABASES_MPI/
MODEL_DIR=../../optimize/
SOLVER_DIR=../../solver/

#parameters to interpolate
MODEL_SET="G0"
KERNEL_SET="dGsp dGcp"

# plot event kernels
PLOT_EVT_KERNEL=0
MDTYPE=dtti
KLTYPE=1

# interpolate options
INTP_KL=0 # if =1 interplate models
INTP_LS=0
run_indx=`seq 0 0`