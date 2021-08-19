#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

echo "Starting to compile surface generator ..."
gcc src/* -lnetcdf -lgeos95 -fopenmp -lm -Wall -o surface_generator
if [ $? -ne 0 ]
then
echo -e ${RED}Surface generator compilation failed.$NC
else
echo "Surface generator compiled successfully."
fi
