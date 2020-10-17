#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

echo "Starting orography generation ..."
mpirun -np $number_of_cpus ./orography_generator $oro_id
# valgrind ./orography_generator $oro_id
if [ $? -ne 0 ]
then
echo -e ${RED}Orography file creation failed.$NC
else
echo "Orography file for oro_id = $oro_id created sucessfully."
fi
