Full Waveform Adjoint Tomography Package for Cube2sph
====================

# specfem3d-cube2sph
download [here](https://github.com/nqdu/specfem3d-cube2sph/tree/devel)

# Source and Stations
1. create a directory `src_rec`, create files called `src_rec/sources.dat.*` (such as `src_rec/sources.dat.noise`), the format is like:
```bash 
NAME evla evlo evdp evbur
```
2. create `src_rec/STATIONS_${NAME}`, the `NAME` should match the first line in `soruces.dat`.

3. `FORCESOLUTION_$NAME` in `src_rec`
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

4. `CMTSOLUTION_$NAME` in `src_rec`


5. data should be in `fwat_data/$NAME`

6. change parameters in `INSTALL` and run it
```bash
./INSTALL slurm INSTALL_DIR
```
then this package will be installed on `INSTALL_DIR`

# python packages required
see `requirements.txt`

