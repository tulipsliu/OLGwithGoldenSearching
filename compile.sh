#!/bin/bash
#Modify the Overlapping Generation Models's Fortran code.
#Author:   Daniel Tulpen Liu
#Date:     July 17, 2024.   Wednesday
#


source /opt/intel/oneapi/setvars.sh --force

ifort -m64 -free -g -qmkl=parallel \
	parameters.f90 global.f90 auxiliary.f90 model_solve.f90 main.f90 tauchen.f90 -o output




