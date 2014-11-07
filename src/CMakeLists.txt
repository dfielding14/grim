# To compiler the code, make a seperate build folder and then do the following:
#
# cd /location/to/build
# cmake /location/to/grim/src
# make
#
# If a different compiler than the one cmake detects is needed, then the
# recommended way is to do the following:
#
# cmake -D CMAKE_C_COMPILER=mpicc -D CMAKE_CXX_COMPILER=mpic++ /location/to/grim/src 

cmake_minimum_required(VERSION 2.8)

project(grim)


# ------------------------------USER OPTIONS----------------------------------#
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3")

# Set custom install folders here
set(PETSC_DIR "/home/mc/Downloads/petsc_optimized")

# Double or float precision?
set(REAL "double")

# Use as ARRAY_ARGS ARRAY_SIZE while declaring arrays in functions arguments
# Must use -std=c99 if using the "restrict" keyword while compiling 
set(ARRAY_ARGS "const restrict static")

# Choose the problem to be solved here:
#
# 1) linear_modes_1d
# 2) komissarov_fast_shock_1d
# 3) komissarov_slow_shock_1d
# 4) komissarov_switch_off_1d
# 5) komissarov_switch_on_1d
# 6) komissarov_shock_tube_1_1d
# 7) komissarov_shock_tube_2_1d
# 8) komissarov_collision_1d
# 9) magnetized_field_loop_advection_2d
# 10) orzag_tang_2d
# 11) komissarov_magnetized_cylindrical_explosion_2d
# 12) equilibrium_torus_2d
# 13) torus_2d
# 14) mti_test_2d
set(PROBLEM "magnetized_field_loop_advection_2d")
# ------------------------------END OF USER OPTIONS---------------------------#

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Need to find all packages before starting compilation, so that all headers
# needed by the code are included.

# Should the external libraries needed by the internal libraries go into their
# own CMakeLists in their own folders?

# TODO: Only find packages depending on what options are set.

# System libraries
find_library(MATH_LIBRARIES m REQUIRED)

# External packages
find_package(PETSc REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(${PETSC_INCLUDES} ${GSL_INCLUDE_DIRS})

add_subdirectory(geometry)
add_subdirectory(physics)
add_subdirectory(reconstruct)
add_subdirectory(riemannsolver)
add_subdirectory(boundary)
add_subdirectory(gridzone)
add_subdirectory(timestepper)

set(PROBLEM_DIR ${CMAKE_SOURCE_DIR}/problem/${PROBLEM})
include(${PROBLEM_DIR}/CMakeLists.txt)
configure_file(${PROBLEM_DIR}/config.h.in ${CMAKE_SOURCE_DIR}/inputs.h)
add_subdirectory(problem)

add_executable(grim grim.c grim.h)

target_link_libraries(grim geometry physics reconstruct boundary gridzone
                      timestepper
                      ${MATH_LIBRARIES}
                      ${PETSC_LIBRARIES})

message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("C Compiler       : " ${CMAKE_C_COMPILER})
message("C_FLAGS          : " ${CMAKE_C_FLAGS})
message("PETSc dir        : " ${PETSC_DIR})
message("")