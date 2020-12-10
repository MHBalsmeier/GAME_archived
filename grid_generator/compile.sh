#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

echo "Starting to compile grid generator ..."
mpicc src/* -O2 -fopenmp -lnetcdf -lm -lgeos95 -Wall -o grid_generator
if [ $? -ne 0 ]
then
echo -e ${RED}Grid generator compilation failed.$NC
else
echo "Grid generator compiled successfully."
fi
