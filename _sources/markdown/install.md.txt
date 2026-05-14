# Introduction

**FWAT-cube2sph** (**F**ull **W**aveform **A**djoint **T**omography) is a package designed for full waveform inversion based on the **specfem3d-cube2sph** solver. 

# Installation 

1. Create Python environment 

we recommend the users to install a new python environment for this package by using `conda` or `mamba`:
```bash
conda create -n fwat python=3.10
conda activate fwat
```
or `venv`:
```bash 
python -m venv fwat 
source fwat/bin/activate
```

2. Install to target directory
```bash
./INSTALL INSTALL_DIR
```
where `INSTALL_DIR` is the absolute path to install your package.
