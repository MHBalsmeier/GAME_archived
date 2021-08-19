#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

echo "***** SURFACE FILE CREATION *****"
if [ $valgrind_check -eq 0 ]
then
mpirun -np $number_of_cpus ./surface_generator $oro_id $oro_rescale_factor
else
if [ $valgrind_check -eq 1 ]
then
valgrind ./surface_generator $oro_id $oro_rescale_factor
fi
fi
if [ $? -ne 0 ]
then
echo -e ${RED}Surface file creation failed.$NC
else
echo "Surface file created sucessfully."
fi
