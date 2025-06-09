Full Waveform Adjoint Tomography Package for Cube2sph
====================

# Source and Stations
1. create a directory `src_rec`, create files called `src_rec/sources.dat.*` (such as `src_rec/sources.dat.noise`), the format is like:
```bash 
NAME evla evlo evdp evbur
```
2. create `src_rec/STATIONS_$NAME`, the `NAME` should match the first line in `soruces.dat`.

3. `FORCESOLUTION/CMTSOLUTION` in `src_rec`

4. data should be in `fwat_data/$NAME`

5. change parameters in `INSTALL` and run it
```bash
./INSTALL INSTALL_DIR
```
then this package will be installed on `INSTAL_DIR`

# python packages required
see `requirements.txt`

