#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

echo "***** SURFACE FILE CREATION *****"

if [ ! -f real/B${res_id}_is_land.nc ]
then
echo "Creating land sea mask ..."
python3 .py/is_land.py $res_id
echo "Land sea mask created."
fi

if [ $valgrind_check -eq 0 ]
then
./surface_generator $oro_id $no_of_cells_for_gliding_avrg
else
if [ $valgrind_check -eq 1 ]
then
valgrind ./surface_generator $oro_id $no_of_cells_for_gliding_avrg
fi
fi
if [ $? -ne 0 ]
then
echo -e ${RED}Surface file creation failed.$NC
else
echo "Surface file created sucessfully."
fi
