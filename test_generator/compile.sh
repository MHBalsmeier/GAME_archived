# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

echo "Starting to compile test generator."
mpicc src/* -leccodes -lnetcdf -lm -lgeos95 -fopenmp -latmostracers -Wl,-rpath=/lib -Wall -o test_generator
if [ $? -ne 0 ]
then
echo -e ${RED}Test generator compilation failed.$NC
else
echo "Test generator compiled sucessfully."
fi
