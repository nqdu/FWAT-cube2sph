# Visualization for Model

This document outlines the usage of visualization scripts for plotting slices of the resulting model. These scripts are located in the plots/plot_model directory. To use them, please ensure that [`GMT`](https://www.generic-mapping-tools.org/) version 6.0 or higher is installed.

---

## Parameter File

This section outlines the parameters used for visualizing slices and event kernels in the `cube2sph` utility. The template is as below:
```bash
#!/bin/bash
specfem_dir="${HOME}/specfem3d-cube2sph/"
cube2sph_dir=$specfem_dir/utils/cube2sph/
NPROC=80

# horizontal slice
NSLICE_HORIZ=1
DEPTH_H=(225)  # add more if you want 
LON0_H=114; LON1_H=132;
LAT0_H=41; LAT1_H=46

# vertical slice
NSLICE_VERTI=1
MAX_DEP=400  # in km

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
```

### Specified Directories
- **`specfem_dir`**: The directory where `specfem3d-cube2sph` is installed.
- **`cube2sph_dir`**: Directory containing the `cube2sph` utility.

### Process Configuration
- **`NPROC`**: Total number of MPI processes, it should be the same as the `NPROC` in `Par_file`.

### Horizontal Slice Parameters
- **`NSLICE_HORIZ`**: Number of horizontal slices to be visualized.
- **`DEPTH_H`**: Depths for horizontal slices (in kilometers).
- **`LON0_H`**: Starting longitude for horizontal slices.
- **`LON1_H`**: Ending longitude for horizontal slices.
- **`LAT0_H`**: Starting latitude for horizontal slices.
- **`LAT1_H`**: Ending latitude for horizontal slices.

### Vertical Slice Parameters
- **`NSLICE_VERTI`**: Number of vertical slices to be visualized.
- **`MAX_DEP`**: Maximum depth for vertical slices (in kilometers).

### Start and End Points for Vertical Slices
- **`LON0_V`**: Starting longitude for vertical slices.
- **`LAT0_V`**: Starting latitude for vertical slices.
- **`LON1_V`**: Ending longitude for vertical slices.
- **`LAT1_V`**: Ending latitude for vertical slices.

### File Directories
- **`DATABASE_DIR`**: Path to the database directory.
- **`MODEL_DIR`**: Path to the model directory.
- **`SOLVER_DIR`**: Path to the solver directory.

### Interpolation Parameters
- **`MODEL_SET`**: Model set to use for visualization.
- **`KERNEL_SET`**: Kernel parameters to interpolate for visualization.

### Additional Visualization Options
- **`PLOT_EVT_KERNEL`**: Flag to indicate whether to plot event kernels (0 = no, 1 = yes).
- **`MDTYPE`**: Model type used for visualization.
- **`KLTYPE`**: Kernel type definition for visualization.

### Interpolation Options
- **`INTP_KL`**: Interpolate models flag (0 = no, 1 = yes).
- **`INTP_LS`**: Interpolate line sources flag.

### Run Index
- **`run_indx`**: Sequence of model indices (0 = `M00`) for visualization tasks.

## Configuration File
Before running any scripts, ensure that you modify your `module_env` and `module_env_gmt` files to load the correct packages and environment variables. For `module_env`, please enable all modules used to compile `cube2sph`, including `gcc`, `gfortran`, `mpi`, and `netcdf-serial`. In `module_env_gmt`, load all packages that are dependencies for **GMT**.


## Step 0: Run MPI Interpolation
This script, `step0_get_gll_from_sem_mpi.sh`, requires MPI to interpolate the slices defined in `parameters.sh`. Run it on your cluster.

After completing the MPI interpolation, you can run `step0_get_mask.sh`. This script will mask (or not, depending on your input variables) the oceans and lakes with `NaN` values and generate a `.grd` file for further use.

## Step 1: Interpolate Model

Run `step1_interp_model.sh`. This script will loop through all `MODEL_SET` and `KERNEL_SET` values, as well as all slices and `run_indx` defined in `parameters.sh`. It will generate text files in the `profiles/` directory.

## Step 2: Generate GRD File

Run `step2_get_grd`. This script will convert all text files generated in **Step 1** into **GMT** `.grd` files for further use.

## Step 3: Plotting

Run `step3_plot_abs.sh` to plot all models and kernels with their true values. If you prefer to visualize the relative values, run `step3_plot_dm.sh`; this will produce figures showing the relative change for all `run_indx` except for `0`. 

Please note that **Step 2** must be executed for model `0` to generate a `.grd` file necessary for computing relative changes.