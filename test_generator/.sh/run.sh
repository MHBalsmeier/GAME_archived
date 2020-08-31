#!/bin/bash

# This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

# test_ids:
# 0:	standard atmosphere without orography
# 1:	standard atmosphere with Gaussian mountain
# 2:	JW test, dry, balanced
# 3:	JW test, dry, perturbed
# 4:	JW test, moist, balanced
# 5:	JW test, moist, perturbed

echo "Starting the test state generation ..."
mpirun -np $number_of_cpus ./test_generator $test_id
# valgrind ./test_generator $test_id
if [ $? -ne 0 ]
then
echo -e ${RED}Test state file creation failed.$NC
else
echo "Test state file for test_id=$test_id created sucessfully."
fi
