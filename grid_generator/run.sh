# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

oro_id=1
echo "Starting grid file creation."
echo "Setup: oro_id = $oro_id."
mpirun -np 1 ./grid_generator $oro_id
if [ $? -ne 0 ]
then
echo -e ${RED}Grid file creation failed.$NC
else
echo "Grid file created sucessfully."
fi
