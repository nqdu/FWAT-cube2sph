# Full Waveform Inversion Tutorial

## FWAT Parameter files
The `fwat_params` directory contains two parameter files:

- **`fwat.yaml`** – Defines the measurement parameters, simulation settings, and optimization options. This file is not automatically refreshed.
- **`lbfgs.yaml`** – Stores the FWI model parameters and is automatically updated after each iteration.

Backup versions of both files are available as `fwat.yaml.org`, and `lbfgs.yaml.org`.

Here is a template of `fwat.yaml`:
```yaml
# FWAT Package parameters
simulation:
  DUMP_WAVEFIELDS: True

# Measurements block, for computing adjoint source
measure:
  # tele seismic
  tele:
    COMPS: ['Z','R'] # components used 
    CH_CODE: BX   # CH_CODE
    FILTER_BANDS: 
      - [5.,50.]
    TIME_WINDOW: [5.,45.] # before and after first arrival
    VERBOSE_MODE: True
    ADJSRC_TYPE: 2 # 2,cross-conv

  noise: # multichannel noise 
    CC_COMPS: ['ZZ']  # {SOURCE-COMP}{RECEIVER-COMP}
    CH_CODE: BX   # CH_CODE
    FILTER_BANDS:  # filter bands used, in s, [T_min,Tmax]
      - [20.,40.]
      - [15.,30.]
      - [10.,20.]
      - [5.,15.]
    GROUPVEL_WIN:   # Time window, determined by group velocity, km/s, [v_min,vmax]
      - [2.0,5.]
      - [2.0,5.]
      - [2.0,5.]
      - [2.0,5.]
    SNR_THRESHOLD: [0.,0.,0.,0.] # exclude data when SNR < SNR_THRESHOLD in a each band 
    USE_EGF: True   # if False, the input data is Cross-correlation, a negative derivative will be applied
    ADJ_SRC_NORM: False  # if true, the adjoint source will be normalized
    USE_NEAR_OFFSET: True # if FALSE, reset tstart
    VERBOSE_MODE: True
    ADJSRC_TYPE: 5 # 5/7/exp_phase/cc_time

  # sks 
  sks:
    COMPS: ['R','T'] # components used 
    CH_CODE: BX   # CH_CODE
    FILTER_BANDS: 
      - [5.,50.]
    TIME_WINDOW: [5.,45.] # before and after first arrival
    VERBOSE_MODE: True
    ADJSRC_TYPE: SI #  SI (splitting_intensity),cross-conv
  
  # receiver function
  rf:
    CH_CODE: BX   # CH_CODE
    FILTER_BANDS: 
      - [5.,50.]
    GAUSS_F0:
      - [1.0]
    MINDERR: 0.001
    MAXIT: 150
    TIME_WINDOW: [5.,25.] # time window, before/after t= 0
    TSHIFT: 5.
    VERBOSE_MODE: True
    ADJSRC_TYPE: 2   # only L2 norm

# optimization
optimize:
  SMOOTHING: [16000.,8000.]  # gaussian smoothing in horizontal/vertical direction, in m
  OPT_METHOD: LBFGS  # GD/ LBFGS
  PRECOND_TYPE: z_precond # default / z_precond /z2_precond
  MAX_PER: 0.02 # maximum relative perturbation

  # model/kernel type
  # iso configuration
  #     0: kappa,mu,rho
  #     1: vp,vs,rho
  #     2: vp/vs,vs,rho
  # dtti configuration:
  #     0: c11-c66,rho
  #     1: vp,vs,rho,gcp,gsp
  #     2. vph,vpv,vsh,vsv,rho,eta,gcp,gsp
  MODEL_TYPE: iso 
  KERNEL_SET: 2
  MASK_VARS: [] # set all the gradients related to index in it to 0
```

several configuration blocks:

### Simulation Block

This block contains parameters for forward and adjoint simulations.

- **`DUMP_WAVEFIELDS`** – If set to `true`, enables the `SUBSAMPLE_FORWARD_WAVEFIELD` option in `Par_file`.  
  In this mode, the full wavefield from the forward simulation is saved and later read back during the adjoint simulation.  

To control how many time steps the wavefield is dumped, modify the following parameters in `Par_file`:
```fortran
KERNEL_SPP = 8
KERNEL_T0  = 5.0
```
- **`KERNEL_SPP`** – Number of sampling points per period.  
- **`KERNEL_T0`** – Minimum period used in your simulation (can be found in `output_generate_databases.txt`).

(measurement-block)=
### Measurement block
This block defines parameters used for measurements, including computing misfits, generating adjoint sources, and applying seismogram rotations.  It currently supports the following four FWI workflows:  

-   Teleseismic waveform inversion  
-   SKS SI-splitting FWI  
-   Ambient noise (single-channel and multi-channel) FWI  
-   Receiver function FWI
Each sub-block specifies settings for a particular measurement method. 

**Note:**  
The `ADJSRC_TYPE` parameter now supports both the [measure_adj](https://github.com/SPECFEM/specfem3d/tree/master/utils/ADJOINT_TOMOGRAPHY_TOOLS/measure_adj) input arguments (1–8) as well as text-based adjoint source types.


#### Teleseismic Waveform Inversion (`tele`)
- **`COMPS`** – List of components used.  
  Example: `['Z','R']`
- **`CH_CODE`** – Channel code.  
  Example: `BX`
- **`FILTER_BANDS`** – Frequency filter bands in seconds, `[T_min, T_max]`.  
  Example: `[[5., 50.]]`
- **`TIME_WINDOW`** – Time window (in seconds) before and after the first arrival.  
  Example: `[5., 45.]`
- **`VERBOSE_MODE`** – If `true`, enables verbose output. Currently it will do nothing.
- **`ADJSRC_TYPE`** – Measurement type code (`2` = L2 norm, or  `cross-conv` = cross convolution). Refer to [measure_adj](https://github.com/SPECFEM/specfem3d/tree/master/utils/ADJOINT_TOMOGRAPHY_TOOLS/measure_adj) for further information.

#### Multi-channel Ambient Noise (`noise`)
- **`CC_COMPS`** – Cross-correlation components in `{SOURCE-COMP}{RECEIVER-COMP}` format. The station files used can be refered to {source section}`source-and-stations`. Example: `['ZZ','TT']`
- **`CH_CODE`** – Channel code.  
  Example: `BX`
- **`FILTER_BANDS`** – List of filter bands in seconds, `[T_min, T_max]`.  
  Example:
`[[20., 40.],
[15., 30.],
[10., 20.],
[5., 15.]]
`
- **`GROUPVEL_WIN`** – Time window determined by group velocity in km/s, `[v_min, v_max]`.  
Example:
`
[[2.5, 4.5],
[2.5, 4.5],
[2.5, 4.5],
[2.5, 4.5]]
`
- **`SNR_THRESHOLD`** - SNR threshold for each frequency band
- **`USE_EGF`** – If `false`, the input data is cross-correlation; a negative derivative will be applied.
- **`ADJ_SRC_NORM`** – If `true`, normalizes the adjoint source.
- **`USE_NEAR_OFFSET`** – If `false`, resets `tstart`.
- **`VERBOSE_MODE`** – If `true`, enables verbose output.
- **`ADJSRC_TYPE`** – Measurement type code (`5` = cross-correlation, or `7` = multitaper, `exp_phase` = exponentiated phase), `cc_time` = 
cross-correlation time misfit

#### SKS SI-Splitting FWI (`sks`)
- **`COMPS`** – List of components used.  
Example: `['R','T']`
- **`CH_CODE`** – Channel code.  
Example: `BX`
- **`FILTER_BANDS`** – Frequency filter bands in seconds, `[T_min, T_max]`.  
Example: `[[5., 50.]]`
- **`TIME_WINDOW`** – Time window (in seconds) before and after the first arrival.  
Example: `[5., 45.]`
- **`VERBOSE_MODE`** – If `true`, enables verbose output.
- **`ADJSRC_TYPE`** – Measurement type code (`SI` = splitting intensity, `cross-conv`).


#### Receiver Functions (`rf`)
- **`CH_CODE`** – Channel code.  
Example: `BX`
- **`FILTER_BANDS`** – Frequency filter bands in seconds, `[T_min, T_max]`.  
Example: `[[5., 50.]]`
- **`GAUSS_F0`** – Gaussian filter width.  
Example: `[[1.0]]`
- **`MINDERR`** – Minimum deconvolution error tolerance.
- **`MAXIT`** – Maximum number of iterations for deconvolution.
- **`TIME_WINDOW`** – Time window (in seconds) before and after `t = 0`.  
Example: `[5., 25.]`
- **`TSHIFT`** – Time shift applied (in seconds).
- **`VERBOSE_MODE`** – If `true`, enables verbose output.
- **`ADJSRC_TYPE`** – Measurement type code (`2` = L2 norm).

### Optimization block
The `optimize` block defines parameters for the optimization process in FWI.

#### General Optimization Settings
- **`SMOOTHING`** – Gaussian smoothing lengths (in meters) for the horizontal and vertical directions.  
  Example: `[16000., 8000.]`
- **`OPT_METHOD`** – Optimization method.  
  Options:  
  - `GD` – Gradient Descent  
  - `LBFGS` – Limited-memory Broyden–Fletcher–Goldfarb–Shanno (default: `LBFGS`)
- **`PRECOND_TYPE`** – Preconditioning method.  
  Options:  
  - `default` – No special preconditioning  
  - `z_precond` – Depth-based preconditioning  
  - `z2_precond` – `z^2` depth-based preconditioning
- **`MAX_PER`** – Maximum relative perturbation allowed during model updates.  
  Example: `0.02` (2% maximum change per iteration)


#### Model and Kernel Settings
- **`MODEL_TYPE`** – Model type.  
  Possible configurations:  

  **Isotropic (`iso`):**
  - `0`: `kappa`, `mu`, `rho`  
  - `1`: `vp`, `vs`, `rho`  
  - `2`: `vp/vs (kappa)`, `vs`, `rho`  

  **Tilted Transversely Isotropic (`dtti`):**
  - `0`: `c11`–`c66`, `rho`  
  - `1`: `vp`, `vs`, `rho`, `gcp`, `gsp`  
  - `2`: `vph`, `vpv`, `vsh`, `vsv`, `rho`, `eta`, `gcp`, `gsp`  

- **`KERNEL_SET`** – Kernel set index.  
  Example: `2`
- **`MASK_VARS`** – List of variable indices to mask in the gradient (set gradients for these indices to `0`).  
  Example: `[]` (no masking)

### LBFGS file
here is a template:
```yaml
# Preconditioned L-BFGS Parameters
MAXITER:  10000   # max iterations
MSTORE: 10         # maximum number of stored pairs
CONV: 1.0e-8

# iteration flag
iter: 0
iter_start: 0
iter_ls: 0
first_ls: True

# FLAG: str one of['INIT','GRAD','PREC','CONV','NSTE','FAIL']
#     = 'INIT' must be used for first iteration
#     = 'GRAD' the user must compute the cost and (preconditioned) gradient at current point x.  
#     = 'PREC' the user must multiply the vector self.q_plb by its preconditioner.  
#     = 'CONV' a minimizer has been found.  
#     = 'NSTE' a new step is performed.    
#     = 'FAIL' the linesearch has failed. 
flag: 'INIT'   

# line search
M1: 1.0e-4 # Wolfe conditions parameter 1 (Nocedal value)
M2: 0.9  # Wolfe conditions parameter 2 (Nocedal value)
FACTOR: 10 # Bracketting parameter (Gilbert value
MAXLS: 100
alpha_L: 0
alpha_R: 0
alpha: -1.
alpha_init: -1.

# debug 
PRINT:  True 
DEBUG: False
```
This block defines settings for the preconditioned L-BFGS optimization algorithm.

At current stage, only part of these parameters are enabled:

#### General Settings
- **`MSTORE`** – Maximum number of stored vector pairs (used in the L-BFGS memory).  
  Example: `10`
- **`iter`** – Current iteration counter.  
- **`iter_start`** – Starting iteration index. Only the memory between `iter_start` and `iter` will be accessed.
- **`iter_ls`** – Line search iteration counter.  

#### Iteration Flag (`flag`)
String flag that indicates the current optimization state.  
Possible values:
- **`INIT`** – First iteration (initialization step).  
- **`GRAD`** – Compute the cost and (preconditioned) gradient at the current point `x`.  
- **`LS`** – line search stage

#### Line Search Parameters
- **`M1`** – Wolfe condition parameter 1 (Nocedal’s value).  
  Example: `1.0e-4`
- **`M2`** – Wolfe condition parameter 2 (Nocedal’s value).  
  Example: `0.9`
- **`FACTOR`** – Bracketing parameter (Gilbert’s value).  
  Example: `10`
- **`alpha_L`** – Left bound for step size.
- **`alpha_R`** – Right bound for step size.
- **`alpha`** – Current step size. For first iteration it wil be `-1.`, then it will be automatically tuned.

(source-and-stations)=
## Source and Stations

This section describes the setup process for source and station files.

### Create the `src_rec` Directory and Source Files
Create a directory named `src_rec` and add source definition files named `src_rec/sources.dat.*`  
(for example, `src_rec/sources.dat.noise`).  

The file format is as follows:
```bash
NAME evla evlo evdp evbur
```
where:
- `NAME` – Source name identifier.
- `evla`, `evlo` – Event latitude and longitude.
- `evdp` – Event depth.
- `evbur` – Event burial depth.

### Create the Station File
Create a file named `src_rec/STATIONS_${NAME}_globe`, where `NAME` matches the **first column** in `sources.dat`.  

Convert it to spherical coordinates using:
```bash
utils/cube2sph/bin/write_station_file
```
This will produce:
```
src_rec/STATIONS_${NAME}, rot_${NAME}
```

### Create the Force Solution File
Create `src_rec/FORCESOLUTION_${NAME}_globe` with the following format:
```
FORCE 001 
time shift:     0.0000
f0:             1.0
latorUTM:       34.5
longorUTM:      125.4
depth:          0.0000
source time function:            0
factor force source:             1.e15
component dir vect source E:     0.e0
component dir vect source N:     0.e0
component dir vect source Z_UP:  1.e0
```
Convert it to:
```
src_rec/FORCESOLUTION_${NAME}
```
using:
```bash
utils/cube2sph/bin/write_force_solution_file
```

**Note:**  
For multi-channel noise simulation, provide the following files:
- `FORCESOLUTION_${NAME}_Z`
- `FORCESOLUTION_${NAME}_N`
- `FORCESOLUTION_${NAME}_E`
- `STATIONS_${NAME}_R` / `STATIONS_${NAME}_T`
- `STATIONS_${NAME}_Z`
The stations in each `STATION_${NAME}_[RTZ]` can be different. You should merge all `rot_${NAME}_{RTZ}` together by 
```bash
cat rot_${NAME}_{RTZ} |sort -n |uniq > rot_${NAME}
```

### Create the CMT Solution File
Create `src_rec/CMTSOLUTION_${NAME}_globe` and convert it to:
```
src_rec/CMTSOLUTION_${NAME}
```
using:
```bash
utils/cube2sph/bin/write_cmt_solution
```

### Place the observation Data
Store the data in:
```
fwat_data/$NAME
```
For multi-channel noise data, it should be placed as:
```
fwat_data/$NAME_[RTZ]
```

## Cluster Parameters

All of the following files are stored in `INSTALL_DIR` during installation.

### `module_env`
Contains the environment module commands required to load dependencies on the cluster.


### `parameters.sh`

- **`SEM_PATH`** – Path to the solver binary directory.  
  Example:  
  ```bash
  SEM_PATH=~/specfem3d-cube2sph
  ```

- **`MPIRUN`** – Command for running MPI jobs.  
  Example:  
  ```bash
  MPIRUN=mpirun
  ```

- **`PLATFORM`** – Execution platform.  
  Options:  
  - `local` – Run locally.  
  - `slurm` – Run on a SLURM-based cluster.  

- **`SIMU_TYPES`** – List of simulation types to run.  
  Example:  
  ```bash
  SIMU_TYPES=("noise")
  ```

- **`SIMU_TYPES_USER_WEIGHT`** – User-defined weights for each simulation type.  
  Example:  
  ```bash
  SIMU_TYPES_USER_WEIGHT=(1.)
  ```

- **`NJOBS_PER_JOBARRAY`** – Number of jobs per job array for each simulation type.  
  Example:  
  ```bash
  NJOBS_PER_JOBARRAY=(1)
  ```

### `run_fwi.sh` and `run_forward.sh`
These scripts implement several functions that add SLURM job headers, making them executable on a cluster.  
The functions include:
- `SHELL_HEADER_SEM`
- `SHELL_HEADER_POST`
- `SHELL_HEADER_WOLFE`
- `WAIT_FINISH`

**Note:**  
Edit these functions to match your working environment, for example, by adding GPU-specific flags, account information, or other SLURM parameters.

## Data Preparation
### Teleseismic Waveform Data
All observed waveforms must be provided in **SAC** format, including all necessary headers (e.g., distance, source and receiver locations, with `lcalda = 1`). Both **R** and **Z** components should be prepared. The reference time (`t = 0`) should be aligned with the theoretical travel times.

under construction

## Checklist

Before running a forward or adjoint simulation, ensure the following items are prepared:

- **`DATA`** – Directory containing input data. You can copy the `DATA` directory used in mesh generation.  
- **`DATA/Par_file.{simu_type}`** – For each simulation type (see [Measurement](#measurement-block)), provide a corresponding `Par_file`.  
- **`DATA/axisem/SOURCE_TAG`** – Background wavefield file, which can be generated using [AxiSEMLib](https://github.com/nqdu/AxiSEMLib/tree/main).  
- **`OUTPUT_FILES`** – You can copy these from the mesh generation output.  
- **`DATABASES_MPI`** – Required for forward simulation.  
- **`./optimize/MODEL_M00`** – Directory containing your initial model (files in `DATABASES_MPI`) for FWI.
- **`fwat_data`** - Data directory. All files must be provided in SAC format with all required headers. Seismograms for each event should be saved in a {NAME} directory, or in {NAME}_{RTZ} for multi-channel noise. For teleseismic data, the time t = 0 must correspond to the direct arrival time.

## User-Defined Parameter Set
under construction

## Visualization 
under construction