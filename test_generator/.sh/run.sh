#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

# test_ids:
# 0:	standard atmosphere without orography
# 1:	standard atmosphere with Gaussian mountain
# 2:	JW test, dry, balanced
# 3:	JW test, dry, perturbed
# 4:	JW test, moist, balanced
# 5:	JW test, moist, perturbed

echo "***** TEST FILE CREATION *****"
echo "Copyright (C) 2020 The GAME development team."
if [ $valgrind_check -eq 0 ]
then
mpirun -np $number_of_cpus ./test_generator $test_id
fi
if [ $valgrind_check -eq 1 ]
then
valgrind ./test_generator $test_id
fi
if [ $? -ne 0 ]
then
echo -e ${RED}Test state file creation failed.$NC
else
echo "Test state file for test_id = $test_id created sucessfully."
fi
