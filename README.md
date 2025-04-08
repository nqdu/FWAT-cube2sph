FWAT Package Usage
====================

# Source and Stations
1. create a directory `src_rec`, create files called `src_rec/sources.dat.*` (such as `src_rec/sources.dat.noise`), the format is like:
```bash 
NAME evla evlo evdp evbur
```
2. create `src_rec/STATIONS_$NAME`, the `NAME` should match the first line in `soruces.dat`.

3. `FORCESOLUTION/CMTSOLUTION` in `src_rec`

4. data should be in `fwat_data/$NAME`

# Compile `measure_adj`
source `module_env` and then go to `measure_adj`, run :
```bash 
make -f Makefile_ifort
```

then go to `fwatlib`, run `compile_all.sh`

# python packages required
see `requirements.txt`

