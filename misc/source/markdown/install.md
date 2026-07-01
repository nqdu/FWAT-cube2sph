# Introduction

**FWAT-cube2sph** (**F**ull **W**aveform **A**djoint **T**omography) is a package for full-waveform inversion based on the **specfem3d-cube2sph** solver.

# Installation 

This package has three parts: a spectral-element solver, the cube2sph toolkit, and the full-waveform inversion (FWI) framework.

1. Create Python environment

We recommend creating a dedicated Python environment for this package using `conda` or `mamba`:
```bash
conda create -n fwat python=3.10
conda activate fwat
```
or `venv`:
```bash 
python -m venv fwat 
source fwat/bin/activate
```

2. Compile the spectral-element solver

First, clone the repo to your installation path:
```
git clone https://github.com/nqdu/specfem3d-cube2sph.git
cd specfem3d-cube2sph
```

Then compile it like:
```bash 
mkdir -p build
cd build
cmake .. -DCC=gcc -DMPIFC=mpif90 
make -j8
```
If you want to build the CUDA-accelerated version, please use:
```bash 
cmake .. -DCC=gcc -DMPIFC=mpif90  -DENABLE_CUDA=ON
```
By default, the CUDA architecture is set to `native`. To target a specific [compute capability](https://developer.nvidia.com/cuda-gpus), open `CMakeLists.txt`, find `set_target_properties(cuda PROPERTIES CUDA_ARCHITECTURES native)`, and replace `native` with your target architecture.

The CUDA-Aware MPI technique would facilitate communications. If you want to enable CUDA-aware MPI, please use:
```bash 
cmake .. -DCC=gcc -DMPIFC=mpif90 -DENABLE_CUDA=ON \
        -DENABLE_CUDA_AWARE=ON
```

3. Compile the cube2sph toolkit

Ensure [netcdf-fortran](https://docs.unidata.ucar.edu/netcdf-fortran/current/) is installed on your machine. Then navigate to `utils/cube2sph` (inside the `specfem3d-cube2sph` repo cloned above) and run:

```bash 
mkdir -p build
cd build
cmake .. -DCC=gcc -DMPIFC=mpif90 
make -j8
```

4. Install the full-waveform inversion scripts

Clone the repo and run the installer:
```bash
git clone https://github.com/nqdu/FWAT-cube2sph.git
cd FWAT-cube2sph
./INSTALL INSTALL_DIR
```
where `INSTALL_DIR` is the absolute path where you want to run FWI.
