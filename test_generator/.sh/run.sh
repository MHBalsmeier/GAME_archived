#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/AUN4GFD/game

echo "***** TEST FILE CREATION *****"
if [ $valgrind_check -eq 0 ]
then
mpirun -np $number_of_cpus ./test_generator $test_id $orography_layers $type_of_vertical_grid $toa
else
if [ $valgrind_check -eq 1 ]
then
valgrind ./test_generator $test_id $orography_layers $type_of_vertical_grid $toa
fi
fi
if [ $? -ne 0 ]
then
echo -e ${RED}Test state file creation failed.$NC
else
echo "Test state file for test_id = $test_id created sucessfully."
fi
