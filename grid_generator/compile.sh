#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

echo "Starting to compile grid generator ..."
gcc -O2 src/* ../src/constituents/dictionary.c -fopenmp -lnetcdf -lm -lgeos95 -Wall -o grid_generator
if [ $? -ne 0 ]
then
	echo -e ${RED}Grid generator compilation failed.$NC
else
	echo "Grid generator compiled successfully."
fi
