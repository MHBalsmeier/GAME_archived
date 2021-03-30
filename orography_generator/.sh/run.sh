#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/AUN4GFD/game

echo "***** OROGRAPHY FILE CREATION *****"
if [ $valgrind_check -eq 0 ]
then
mpirun -np $number_of_cpus ./orography_generator $oro_id $rescale_factor
else
if [ $valgrind_check -eq 1 ]
then
valgrind ./orography_generator $oro_id $rescale_factor
fi
fi
if [ $? -ne 0 ]
then
echo -e ${RED}Orography file creation failed.$NC
else
echo "Orography file for oro_id = $oro_id created sucessfully."
fi
