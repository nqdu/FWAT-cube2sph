# Examples

This section provides step-by-step examples for using the **FWAT-cube2sph** package, including mesh generation, forward simulation, and full-waveform inversion (FWI).

## Example: FWI with a Checkerboard Model (NED-based)

This example demonstrates how to set up a complete FWI workflow using a checkerboard velocity model based on the NED (National Elevation Data) dataset. The workflow includes mesh generation, forward simulation to create synthetic data, and FWI inversion.

### Prerequisites

Ensure you have completed the installation steps from the [Installation](install.md) guide. The FWI workflow is typically set up in a dedicated project directory (e.g., `NED/fwi/`):

```bash
NED/fwi/
├── DATA/                          # SPECFEM input data directory
├── src_rec/                       # Source and station definitions
├── fwat_data/                     # Observed/synthetic seismic data
├── optimize/                      # Initial model and optimization outputs
├── config.env                     # Environment configuration
├── fwat_params/
│   ├── fwat.yaml                 # FWI parameters
│   └── lbfgs.yaml                # Optimization parameters
├── run_forward.sh                 # Forward simulation script
└── run_fwi.sh                     # FWI workflow script
```

The `src_rec/` and `fwat_params/` directories are pre-configured in your project directory. Customize them as needed for your inversion.

### Step 1: Mesh Generation

The first step is to create a structured mesh based on your velocity model.

#### 1a. Create the Checkerboard Model

Generate a checkerboard velocity model using the provided utilities:

```bash
cd mesh/create_model
python tomo.py        # Create background tomography model
python add_ckbd_model.py  # Add checkerboard perturbations
```

This produces a `tomography_model.xyz` file that defines the velocity model.

#### 1b. Create Tag-specific Models

For each simulation type (e.g., `1d`, `ckbd`), prepare separate model files and generate corresponding mesh databases:

```bash
cd ../..
mkdir -p optimize/tomo.1d
mkdir -p optimize/tomo.ckbd

# Process 1D model
cp mesh/create_model/tomography_model.xyz.1d mesh/DATA/tomo_files/tomography_model.xyz
cd mesh
bash step1-step2     # Generate mesh and databases
cp -r DATABASES_MPI OUTPUT_FILES ../optimize/tomo.1d/

# Process checkerboard model
cd create_model
cp tomography_model.xyz.ckbd ../DATA/tomo_files/tomography_model.xyz
cd ..
bash step1-step2
cp -r DATABASES_MPI OUTPUT_FILES ../optimize/tomo.ckbd/
cd ..
```

The `step1-step2` script performs:
- Step 1: Generate mesh (`xgenerate_databases`)
- Step 2: Prepare mesh for FWI

After completion:
- `DATABASES_MPI` – Contains the mesh database files
- `OUTPUT_FILES` – Contains mesh and simulation configuration outputs

### Step 2: Customize Source and Station Configuration

Source and station definitions are pre-configured in the `src_rec/` directory. Customize them for your inversion setup.

#### 2a. Source Definitions

Edit `src_rec/sources.dat.noise` to define your sources:

```bash
NOISE1  35.0  120.0  10.0  0.0
```

Format: `NAME  latitude  longitude  depth  burial_depth`

#### 2b. Station Definitions

Create station files in geographic coordinates (e.g., `src_rec/STATIONS_NOISE1_globe`) and convert to cube coordinates:

```bash
utils/cube2sph/bin/write_station_file
```

This generates `src_rec/STATIONS_NOISE1` and `src_rec/rot_NOISE1`.

#### 2c. Source Time Functions

Define force solution files in geographic coordinates (e.g., `src_rec/FORCESOLUTION_NOISE1_Z_globe`) and convert to cube coordinates:

```bash
utils/cube2sph/bin/write_force_solution_file
```

Create files for each component (Z, N, E).

### Step 3: Configure FWI Parameters

Edit `fwat_params/fwat.yaml` to configure the inversion. Refer to the [Full Waveform Inversion Tutorial](fwat.md) for detailed parameter descriptions.

### Step 4: Set Up Cluster Configuration

Edit `config.env` with your system environment:

```bash
# Module loading (if needed)
# module load gcc/9.3
# module load openmpi/4.0

# Path to SPECFEM solver
SEM_PATH=/path/to/specfem3d-cube2sph

# MPI command
MPIRUN=mpirun

# Number of processes for measurement
NPROC_MEASURE=16

# Platform: local, slurm, or pbs
PLATFORM="local"

# Local execution settings
USE_IO_TMPDIR=0
```

### Step 5: Forward Simulation

Generate synthetic data using the checkerboard model (approximately 10 seconds on RTX 4070):

```bash
bash run_forward.sh
```

This step:
1. Runs the forward simulation for each source
2. Produces synthetic seismic data (wavefields, seismograms)
3. Saves data to `fwat_data/NOISE1/`

### Step 6: Run FWI

Once synthetic data is ready, start the inversion (approximately 4 hours with 1 node and 2 processes per node):

```bash
bash run_fwi_big_job.sh 1 2 35
```

The FWI workflow will:
1. Run forward simulations for each source with the current model
2. Compare synthetic data with observations
3. Compute misfit and gradients
4. Update the model using L-BFGS optimization
5. Repeat for each iteration

### Step 7: Monitor and Analyze Results

#### Check Inversion Progress

```bash
# View misfit values
cat optimize/misfit_iter.txt

# Examine model updates
ls optimize/MODEL_M*/
```

#### Plot Results

To visualize waveforms and model updates, use the plotting utilities:

```bash
# Plot misfit traces
cd plots/plot_traces
# Run relevant plotting scripts

# Plot model evolution
cd plots/plot_model
# Run step0-3 to visualize model updates across iterations
```

### Practical Tips

- **Start with coarse resolution**: Begin with a checkerboard model at low frequency bands to test the workflow before running full-resolution inversions.

- **Monitor runtime**: Check the first iteration's runtime to estimate total inversion time. Adjust `NPROC_MEASURE` or number of nodes if jobs are slow.

- **Adjust smoothing parameters**: If the inversion produces oscillatory models, increase `SMOOTHING` in `fwat.yaml`.

- **Data quality**: Ensure input seismic data are in SAC format with correct headers.

- **Restart from checkpoint**: If an inversion is interrupted, simply rerun `bash run_fwi.sh` to resume from the last saved iteration.

## Next Steps

For more details on FWI parameters, data preparation, and optimization algorithms, refer to the [Full Waveform Inversion Tutorial](fwat.md).
