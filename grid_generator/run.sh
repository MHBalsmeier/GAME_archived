# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

oro_id=2
optimize=1
n_iterations=2000
use_scalar_h_coords_file=1
scalar_h_coords_file="nc_files/B5L26T30000_O0_OL17_SCVT.nc"
echo "Starting grid file creation."
echo "Setup: oro_id = $oro_id."
mpirun -np 1 ./grid_generator $oro_id $optimize $n_iterations $use_scalar_h_coords_file $scalar_h_coords_file
if [ $? -ne 0 ]
then
echo -e ${RED}Grid file creation failed.$NC
else
echo "Grid file created sucessfully."
fi
