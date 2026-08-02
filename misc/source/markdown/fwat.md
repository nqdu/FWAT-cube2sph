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

  types: ["noise"] # simulation types: noise/tele/sks/rf
  weights: [1.] # weight for each simu type, used for multiple simu types, e.g., ["noise","tele"], [0.5,0.5]

  # # < 0 NO ADAPTIVE WEIGHTS, >=0 use this iteration misfit to compute weights chi = \sum chi_i * w_i
  # chi = \sum_i L_i(x) / L_i(x0) / |g_i(x0)| * w_i, and weights will be re-written as (w_i / L_i(x0) / |g_i(x0)|)
  iter_wts: 0 

  # Normalization by gradient norm or misfit at iter_wts, if iter_wts >= 0
  norm_type: 'misfit' # gradient or misfit


# Measurements block, for computing adjoint source
measure:
  # tele seismic
  tele:
    COMPS: ['Z','R'] # components used 
    CH_CODE: BX   # CH_CODE
    FILTER_BANDS: 
      - [5.,50.]
    TIME_WINDOW: [5.,45.] # before and after first arrival
    VERBOSE_MODE: False
    ADJSRC_TYPE: 2 # 2 (l2),cross-conv,cc_time_dd

  noise: # multichannel noise 
    CC_COMPS: ['ZZ','TT']  # {SOURCE-COMP}{RECEIVER-COMP}
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
    TSHIFT_MAX: [4.5,4.5,4.5,4.5] # exclude data when ABS(TSHIFT) > TSHIFT_MAX in a each band
    DLNA_MAX: [1.5,1.5,1.5,1.5] # exclude data when ABS(DLNA) > DLNA_MAX in a each band
    CC_MIN: [0.8,0.8,0.8,0.8] # exclude data when CC_COEF < CC_MIN in a each band
    USE_EGF: True   # if False, the input data is Cross-correlation, a negative derivative will be applied
    ADJ_SRC_NORM: False  # if true, the adjoint source will be normalized
    USE_NEAR_OFFSET: True # if FALSE, reset tstart
    VERBOSE_MODE: True
    ADJSRC_TYPE: 5 # 5/7/exp_phase/cc_time/cc_time_dd

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
  PRECOND_TYPE: none # default / z_precond /z2_precond/ z_sqrt_precond/ none
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
  #     3. vp,vs,rho,kappaa,kappab,eta,gcp,gsp  where kappaa = (vph-vpv)/vpv, vp = sqrt((2*vsh^2 + vsv^2)/3)
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

- **`types`** – List of simulation types to run.  
  Options: `noise`, `tele`, `sks`, `rf`.  
  Example: `["noise"]`
- **`weights`** – User defined Weight for each simulation type, used when multiple types are combined.  
  Example: `[0.5, 0.5]` for `["noise", "tele"]`
- **`iter_wts`** – Controls adaptive weighting. If `< 0`, adaptive weights are disabled. If `>= 0`, the misfit at this iteration is used to compute weights as:  
  $\chi = \sum_i L_i(x) \cdot w_i$  
  and the weights are rewritten as $w_i * norm_i$, where $norm_i$ is the normalization factor  
  Example: `0`
- **`norm_type`** – Normalization strategy when adaptive weights are enabled (`iter_wts >= 0`). Controls whether the weight normalization factor is computed from the gradient norm or the misfit value.  
  Options: `gradient`, `misfit`  
  Example: `'misfit'`

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
- **`SNR_THRESHOLD`** – SNR threshold for each frequency band; data are excluded when SNR < `SNR_THRESHOLD`.
- **`TSHIFT_MAX`** – Maximum allowed absolute time shift (in seconds) per frequency band; data are excluded when |TSHIFT| > `TSHIFT_MAX`.
- **`DLNA_MAX`** – Maximum allowed absolute amplitude ratio (dlna) per frequency band; data are excluded when |DLNA| > `DLNA_MAX`.
- **`CC_MIN`** – Minimum cross-correlation coefficient per frequency band; data are excluded when CC_COEF < `CC_MIN`.
- **`USE_EGF`** – If `false`, the input data is cross-correlation; a negative derivative will be applied.
- **`ADJ_SRC_NORM`** – If `true`, normalizes the adjoint source.
- **`USE_NEAR_OFFSET`** – If `false`, resets `tstart`.
- **`VERBOSE_MODE`** – If `true`, enables verbose output.
- **`ADJSRC_TYPE`** – Measurement type code (`5` = cross-correlation, or `7` = multitaper, `exp_phase` = exponentiated phase), `cc_time` = 
cross-correlation time misfit, `cc_time_dd` = double difference cc time misfit.

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
  - `z_sqrt_precond` – Square-root depth-based preconditioning  
  - `none` – No preconditioning
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
  - `3`: `vp`, `vs`, `rho`, `kappaa`,`kappab`, `eta`, `gcp`,`gsp`

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
- `STATIONS_${NAME}_RR/TT/ZZ/RZ` (stations used for `RR,TT,ZZ,RZ`)

The stations in each `STATION_${NAME}_[RR/TT/ZZ/RT]` can be different. You should merge all `rot_${NAME}_{RTZ}` together by using this python file
```python
import numpy as np
import glob 
import sys 
import os 

def read_rot_file(mydict:dict,filename):

    with open(filename,'r') as fio:
        lines = fio.readlines()
        n = 0
        while n < len(lines):
            line = lines[n].strip()
            
            # get keys
            info = line.split()
            key = info[0] + '_' + info[1]

            # get values, 3x3 matrix
            a = np.zeros((3,3),dtype=np.float64)
            for i in range(3):
                n += 1
                line = lines[n].strip()
                info = line.split()
                for j in range(3):
                    a[i,j] = float(info[j])
            
            n += 1
            # check if key exists
            if key not in mydict:
                mydict[key] = a
    
    return mydict

def main():
    if len(sys.argv) < 2:
        print("Usage: python merge_rot_file.py path ...")
        print("example: python merge_rot_file.py src_rec")
        sys.exit(1)
    
    path = sys.argv[1]

    # get all rot files, start with rot_*_[ZZ/TT/RR]
    input_files = glob.glob(f"{path}/rot_*_*")

    # got all names
    all_names = set()
    for f in input_files:
        parts = os.path.basename(f).split('_')
        all_names.add(parts[1])

    # loop each names
    for name in all_names:
        mydict = dict()
        filenames = glob.glob(f"{path}/rot_{name}_*")
        for f in filenames:
            if f"rot_{name}_" in f:
                mydict = read_rot_file(mydict,f)

                # remove the file after reading
                os.remove(f)
        
        # write merged file
        output_file = f"{path}/rot_{name}"
        with open(output_file,'w') as fio:
            keys=  sorted(mydict.keys())
            for key in keys:
                a,b = key.split('_')
                fio.write(f"{a}\t {b}\n")
                a = mydict[key]
                for i in range(3):
                    fio.write(f"{a[i,0]}\t {a[i,1]}\t {a[i,2]}\n")

if __name__ == "__main__":
    main()

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

### `config.env`
Contains the shell environment setup and runtime options used by the workflow scripts.

- Shell setup lines at the top of this file are executed before the workflow starts. Use them to load modules, source site-specific environment scripts, and export MPI or OpenMP variables required on your system.
- **`SEM_PATH`** – Path to the SPECFEM solver directory. The scripts expect solver binaries such as `xspecfem3D`, `xsmooth_sem_sph_pde`, and `xgenerate_databases` under `$SEM_PATH/bin`.

  Example:
  ```bash
  SEM_PATH=~/software/specfem3d-cube2sph
  ```
- **`MPIRUN`** – Command used to launch MPI jobs.

  Example:
  ```bash
  MPIRUN=mpirun
  ```
- **`NPROC_MEASURE`** – Number of MPI ranks used by the measurement step (`fwat-main measure`). This may differ from the number of solver ranks used by SPECFEM.

  Example:
  ```bash
  NPROC_MEASURE=16
  ```
- **`PLATFORM`** – Execution platform.
  Options:
  - `local` – Run locally.
  - `slurm` – Run on a SLURM-based cluster.
  - `pbs` – Run on a PBS-based cluster.

  Example:
  ```bash
  PLATFORM="slurm"
  ```
- **`USE_IO_TMPDIR`** – If set to `1`, the measurement scripts stage temporary working files in node-local storage when available. Set it to `0` to keep all I/O in the run directory.

  Example:
  ```bash
  USE_IO_TMPDIR=1
  ```
- **`IO_TMPDIR`** – Path to the node-local temporary directory used when `USE_IO_TMPDIR=1`. On SLURM systems this is commonly set to `$SLURM_TMPDIR`.

  Example:
  ```bash
  IO_TMPDIR=$SLURM_TMPDIR
  ```
- **`NJOBS_PER_JOBARRAY`** – Number of jobs per job array for each simulation type. Its length must match `simulation.types` in `fwat.yaml` when running on `slurm` or `pbs`.

  Example:
  ```bash
  NJOBS_PER_JOBARRAY=(1)
  ```

### `run_fwi.sh` and `run_forward.sh`

These scripts implement several functions to add SLURM job headers, enabling execution on a cluster. Both scripts leverage the job array capabilities of SLURM.  

The functions included are:
- **`SHELL_HEADER_SEM`**
- **`SHELL_HEADER_POST`**
- **`SHELL_HEADER_WOLFE`**
- **`WAIT_FINISH`**

**Note:**  
Be sure to modify these functions according to your working environment. This may involve adding GPU-specific flags, account details, or other SLURM parameters.

### `run_fwi_big_job.sh`

This script facilitates an FWI workflow and serves as an alternative to `run_fwi.sh`. It launches a large job that may span multiple nodes and run for several hours. The input parameters specify the number of processes and the number of iterations. Events will be allocated across all processes utilized in this script, and the workflow will continue iterating until it completes the specified iterations or the allotted time expires.  

To submit the job locally, use the following command:
```bash
bash run_fwi_big_job.sh 32 10
```
In this example, the job will launch `32` MPI processes and iterate `10` times.

To run the job on SLURM/PBS, use:
```bash
bash run_fwi_big_job.sh 4 64 10
```
This command will launch a SLURM/PBS job with `4` nodes and `64` MPI processes per node (a total of 256 processes), and will iterate `10` times.

In both cases, the script will launch `NPROCS_IN_TOTAL / NPROCS` simulations (where `NPROCS` is defined in the `Par_file`).

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

If you want to define your own type of elastic model, you can modify `fwat/optimize/model.py` by defining your own `mdtype` and `kltype`. Here are several important notes:

- **Do not change the base model**: For the isotropic case, stick to the parameters `['vp', 'vs', 'rho']`, and for the anisotropic case, the base model should be `[c11-c66, rho]`.
  
- **Define your parameter set**: Set the search direction names in the function `direc_names` by adding a prefix `d` to each of your parameters.
  
- **Implement your own `convert_model`**: This function should convert the base model to your parameter set and vice versa.

- **Implement `convert_to_visual`**: This function will return the visualization model for plotting. It can differ from your parameter set. For example, you may include `gcp` and `gsp` in your parameter set, but for visualization, it is preferable to use `G0` and `phi`.

- **Implement your `user2opt`**: This function maps your user-defined model into the optimization parameter space used by the gradient-based optimizer (e.g., taking `log(vs)` for parameters with physical units, while leaving dimensionless parameters unchanged).

- **Implement `model_update`**: This function guides you in updating your model. For dimensionless parameters, you will directly add the search direction, and for parameters with physical units, you will apply an exponential factor based on the search direction.

- **Implement `convert_kl`**: This function converts the base-model kernels into kernels in the optimization parameter space (consistent with `user2opt`). Do not re-derive these derivatives by hand — for the `dtti` model type they are auto-generated (see the guide below).

### Auto-generating the `dtti` conversions with `auto_kergen.py`

For the `dtti` model type you do **not** hand-write the parameter/kernel conversions. They are generated by `fwat/optimize/auto_kergen.py` into `fwat/optimize/dtti_generated.py`, which `model.py` imports; `convert_model` and `convert_kl` are only thin dispatchers into that generated module.

Everything is derived from a single symbolic stiffness matrix $C_{ij}(\text{params})$ per `kltype`, so the model transform and its kernel chain rule can never drift apart. For each `kltype` the generator emits three things:

- **backward model** (user params → `cijkl`): evaluate $C_{ij}(\text{params})$;
- **forward model** (`cijkl` → user params): the analytic inverse/projection (a modelling choice, e.g. $L = (c_{44}+c_{55})/2$);
- **kernels** (base `cijkl` kernels → optimization space): the chain rule $K_p = \sum_{ij} \dfrac{\partial C_{ij}}{\partial p}\, K_{c_{ij}}$.

`sympy` performs the differentiation and common-subexpression elimination (CSE), and plain `numpy` code is written out.

**Commands**
```bash
cd src/fwat/optimize
# (re)generate dtti_generated.py -- run this after editing auto_kergen.py
python auto_kergen.py --write
# print one kltype kernel body to stdout, just to inspect it
python auto_kergen.py dtti 3
```

**To add a new `kltype`** (all edits are in `auto_kergen.py`):

1. In `_case`, add a branch that names your parameters, builds the 6×6 Voigt stiffness `C0` from them (helper `_c66mat`), and sets `n_log` — the number of leading parameters carried in log space (velocities and density).
2. In `_forward_exprs`, add the matching inverse that recovers your parameters from the base `cijkl` vector.
3. Run `python auto_kergen.py --write` to regenerate `dtti_generated.py`.
4. Add the new `kltype`'s names to `direc_names`, and adjust `user2opt` / `model_update` if your log-space split differs.

**Important**

- Do **not** edit `dtti_generated.py` by hand — your changes are overwritten on the next `--write`.
- Validate with `python tests/test_grad_check.py` (run from the repository root), which finite-difference-checks the analytic kernels against the misfit. A maximum relative error at the ~0.05% level is finite-difference truncation error, not a bug (the dimensionless parameters, which do not depend on the step size, agree to ~1e-6%).