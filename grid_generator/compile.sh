# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

echo "Starting to compile grid generator ..."
mpicc -O3 src/* -lnetcdf -lm -lgeos95 -Wall -o grid_generator
if [ $? -ne 0 ]
then
echo -e ${RED}Grid generator compilation failed.$NC
else
echo "Grid generator compiled successfully."
fi
