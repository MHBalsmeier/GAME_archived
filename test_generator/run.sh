# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

# test_ids:
# 0:	standard atmosphere without orography
# 1:	standard atmosphere with Gaussian mountain
# 2:	JW test, dry, balanced
# 3:	JW test, dry, perturbed
# 4:	JW test, moist, balanced
# 5:	JW test, moist, perturbed


test_id=2
echo "Starting the test state generation ..."
mpirun -np 1 ./test_generator $test_id
# valgrind ./test_generator $test_id
if [ $? -ne 0 ]
then
echo -e ${RED}Test state file creation failed.$NC
else
echo "Test state file created sucessfully."
fi
