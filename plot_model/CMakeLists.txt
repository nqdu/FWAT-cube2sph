#  CMakeLists.txt

cmake_minimum_required(VERSION 3.5.0)

################## compiler options ##################
# compiler
set(CMAKE_Fortran_COMPILER ${FC})

# compiler flags

################## dependencies ##################

# set default values
macro(setDefault var val)
if(NOT DEFINED ${var})
set(${var} ${val})
endif()
endmacro()

# project
project(interp3D  Fortran)

# build type 
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/build/obj)

################# cmake setup #################
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_FORTRAN_FLAGS} -std=f2008")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_FORTRAN_FLAGS} -std08")
endif()
set(CMAKE_VERBOSE_MAKEFILE off)

################# find packages #################
# path of *.cmake files
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/include/cmake)

# mpi
find_package(MPI)
include_directories(${MPI_FORTRAN_INCLUDE_DIRS})
set(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS}")

################# include paths #################
include_directories(${PROJECT_SOURCE_DIR}/src/)

# subs
add_subdirectory(${PROJECT_SOURCE_DIR}/src/interp3D)

# generate_mesh
add_executable(generate_slice ${PROJECT_SOURCE_DIR}/src/interp3D/prog_gen_slices.f90)
add_executable(interp3D ${PROJECT_SOURCE_DIR}/src/interp3D/prog_interp3D.f90)
add_executable(interp3D_serial ${PROJECT_SOURCE_DIR}/src/interp3D/prog_serial.f90)
target_link_libraries(interp3D shared ${MPI_Fortran_LIBRARIES})
target_link_libraries(generate_slice shared ${MPI_Fortran_LIBRARIES})
target_link_libraries(interp3D_serial shared ${MPI_Fortran_LIBRARIES})