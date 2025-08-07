Full Waveform Adjoint Tomography Package for Cube2sph
====================

# specfem3d-cube2sph
download [here](https://github.com/nqdu/specfem3d-cube2sph/tree/devel)

# Source and Stations
1. create a directory `src_rec`, create files called `src_rec/sources.dat.*` (such as `src_rec/sources.dat.noise`), the format is like:
```bash 
NAME evla evlo evdp evbur
```

2. create `src_rec/STATIONS_${NAME}_globe`, the `NAME` should match the first line in `soruces.dat`. Then convert it to `src_rec/STATIONS_${NAME}` using `utils/cube2sph/bin/write_station_file`

3. `FORCESOLUTION_$NAME_globe` in `src_rec`
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
Then convert it to `FORCESOLUTION_$NAME` using `utils/cube2sph/bin/write_force_solution_file`. For multi-channel noise simulation, you should provide `FORCESOLUTION_$NAME_Z`,`FORCESOLUTION_$NAME_N`,`FORCESOLUTION_$NAME_E`, and `STATIONS_${NAME}_R/T`, 
`STATIONS_${NAME}_Z`

4. `CMTSOLUTION_$NAME_globe` in `src_rec`. Then convert it to `CMTSOLUTION_$NAME` using `utils/cube2sph/bin/write_cmt_solution`


5. data should be in `fwat_data/$NAME`

6. change parameters in `INSTALL` and run it
```bash
./INSTALL slurm INSTALL_DIR
```
then the submission scripts will be installed in `INSTALL_DIR`


