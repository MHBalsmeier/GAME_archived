# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

oro_id=2
optimize=1
n_iterations=2000
use_scalar_h_coords_file=1
scalar_h_coords_file="nc_files/B5L26T30000_O0_OL17_SCVT.nc"

# END OF INPUT SECTION

echo "***** GRID FILE CREATION *****"
echo "(C) 2020 The GAME development team."
echo ""
echo "Setup:"
echo "oro_id = $oro_id."
if [ $optimize -eq 1 ]
then
echo "Optimization is switched on. Number of iterations: $n_iterations."
else
echo "Optimization is switched off."
fi
if [ $use_scalar_h_coords_file -eq 1 ]
then
echo "Horizontal coordinates of the generating points (the scalar points in terms of the model) will be read from file $scalar_h_coords_file. This overwrites optimization options."
fi
echo ""
echo "********** Calling the GAME grid generator **********"
echo ""
mpirun -np 1 ./grid_generator $oro_id $optimize $n_iterations $use_scalar_h_coords_file $scalar_h_coords_file
if [ $? -ne 0 ]
then
echo -e ${RED}Grid file creation failed.$NC
else
echo "Grid file created sucessfully."
fi
