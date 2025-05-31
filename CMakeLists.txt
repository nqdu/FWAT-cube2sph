cmake_minimum_required(VERSION 3.12)

# compiler
set(CMAKE_CXX_COMPILER          ${CXX})

# compiler flags
set(CMAKE_Fortran_COMPILER      ${FC}) 
set(SACHOME      ${SACHOME}) 

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_VERBOSE_MAKEFILE ON)

project(measure_adj LANGUAGES CXX Fortran)
file (GLOB SOURCE_FILES "${PROJECT_SOURCE_DIR}/src/*.f90")

#find_package(Python3 COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(libmeas src/main.cpp ${SOURCE_FILES})

#message(${SACHOME})
target_link_directories(libmeas PRIVATE ${SACHOME}/lib)
#link_directories(${SACHOME}/lib)
target_link_libraries(libmeas PRIVATE sacio sac)
#set_target_properties(sacio sac  PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(TARGETS libmeas
  COMPONENT python
  LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/fwatlib/measure/lib
)