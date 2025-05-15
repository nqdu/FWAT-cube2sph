SACHOME=${HOME}/software//sac/
SACLIBDIR=${SACHOME}/lib
LIB="-lsacio -lsac"

set -x 
mpif90 -O3 ma_constants.f90 ma_sub.f90 ma_sub2.f90 ma_variables.f90 ascii_rw.f90 measure_adj_mod.f90 sacio.f90 measure_adj_mpi.f90 -Imod -Jmod -o measure_adj_mpi -L${SACLIBDIR} ${LIB}