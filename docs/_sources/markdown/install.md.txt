# Introduction

**FWAT-cube2sph** (**F**ull **W**aveform **A**djoint **T**omography) is a package designed for full waveform inversion based on the **specfem3d-cube2sph** solver. 

# Installation 

1. **Compilers:** C++/Fortran compilers which support c++17 (tested on `GCC >=7.5`), `cmake >= 3.12`

2. Python environment 

we recommend the users to install a new python environment for this package by using `conda`:
```bash
conda create -n fwat python=3.10
conda activate fwat
```
or `venv`:
```bash 
python -m venv fwat 
source fwat/bin/activate
```

4. Install
change `CXX` and `F90` in `INSTALL`, and run it like
```bash
./INSTALL INSTALL_DIR
```
where `INSTALL_DIR` is the absolute path to install your package.
