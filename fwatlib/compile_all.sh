#!/bin/bash

# # gfortran
# CC="gcc -O3 -fPIC -shared"
# FC="gfortran -O3 -fPIC -shared"
# linkf="-lgfortran -lopenblas"

#ifort
CC="icc -O3 -fPIC -shared"
FC="ifort -O3 -fPIC -shared"
linkf="-lifcore -lmkl_rt  -lpthread -lm -ldl"

#FC="gfortran  -fPIC"
include=`python -m pybind11 --includes`
outname="libpca"`python3-config --extension-suffix`

chmod +x *.sh 

set -x 
cd measure/tele/pca 
$FC -g -c spanlib.f90 -o spanlib.o
$FC -g -c pcalib.f90 -o pcalib.o
$CC -g -c main.cpp -o main.o $include
$CC *.o -o ./$outname $linkf  

rm *.o *.mod
