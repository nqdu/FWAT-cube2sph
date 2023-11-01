file(GLOB_RECURSE SOURCE *.f90)
list(FILTER SOURCE EXCLUDE REGEX "prog")
message(${SOURCE})
add_library(shared ${SOURCE})