#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

echo "Starting to compile orography generator ..."
gcc src/* ../shared/various.c -lnetcdf -lgeos95 -fopenmp -lm -Wall -o orography_generator
if [ $? -ne 0 ]
then
echo -e ${RED}Orography generator compilation failed.$NC
else
echo "Orography generator compiled successfully."
fi
