#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

echo "***** OROGRAPHY FILE CREATION *****"
echo "Copyright (C) 2020 The GAME development team."
if [ $valgrind_check -eq 0 ]
then
mpirun -np $number_of_cpus ./orography_generator $oro_id
fi
if [ $valgrind_check -eq 1 ]
then
valgrind ./orography_generator $oro_id
fi
if [ $? -ne 0 ]
then
echo -e ${RED}Orography file creation failed.$NC
else
echo "Orography file for oro_id = $oro_id created sucessfully."
fi
