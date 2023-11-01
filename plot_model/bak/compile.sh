#/bin/bash
module load intel openmpi

FC="mpif90 -std08 -traceback -warn all "
cd src 
set -e -x

$FC -g -c utm_geo.f90 -o utm_geo.o -O3 
$FC -g -c read_media3D_mpi.f90 -o read_media3D.o -O3 
$FC *.o -o ../bin/read_media 
rm *.o 
cd ..

# module load gcc/8.3.0
# CXX=g++
# FC=gfortran
# CXXFLAGS="-g -c -std=c++14 -O3 -fPIC -shared"
# FFLAGS="-g -c -std=f2008 -fPIC -shared -O3 "
# PYBIND11_INC=$(python -m pybind11 --includes)
# dllname="../libutm"$(python3-config --extension-suffix)

 
# set -x 
# cd src 
# $CXX $CXXFLAGS utm_projection.cpp -o utm_projection.o  $PYBIND11_INC
# $FC $FFLAGS utm_geo.f90 -o utm_geo.o 
# $CXX *.o -o $dllname -lgfortran -fPIC -shared
# rm *.o 
# cd .
