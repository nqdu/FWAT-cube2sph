#/bin/bash
FC="mpiifort -std08 -check nobounds"
cd src 
set -x 

$FC -g -c utm_geo.f90 -o utm_geo.o -O3 
$FC -g -c read_media3D_mpi.f90 -o read_media3D.o -O3 
$FC *.o -o ../bin/read_media 
rm *.o 
cd ..
