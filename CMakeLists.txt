cmake_minimum_required(VERSION 3.5)

project(asteroid2022)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -O2 -std=f2008 -fimplicit-none -Wextra -Wall -pedantic")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone -warn all")
    set(bounds "-check bounds")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

#
# Compile.
#
file(GLOB_RECURSE sources  src/*.f90 src/*.h)
add_executable(prog ${sources})