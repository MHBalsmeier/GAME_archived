#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# verbosity
echo "***** GRID FILE CREATION *****"
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
if [ ! -f $scalar_h_coords_file ]
then
echo "$scalar_h_coords_file does not exist."
echo "Aborting."
fi
echo "Horizontal coordinates of the generating points (the scalar points in terms of the model) will be read from file $scalar_h_coords_file."
fi
echo "number of layers following orography: "$orography_layers
echo "stretching parameter: "$stretching_parameter
echo "model top: "$toa" m"
echo "type of vertical grid: "$type_of_vertical_grid
# end verbosity

echo ""
echo "********** Calling the GAME grid generator **********"
echo ""
if [ $valgrind_check -eq 0 ]
then
./grid_generator $oro_id $optimize $n_iterations $use_scalar_h_coords_file $scalar_h_coords_file $stretching_parameter $orography_layers $type_of_vertical_grid $toa
else
if [ $valgrind_check -eq 1 ]
then
valgrind ./grid_generator $oro_id $optimize $n_iterations $use_scalar_h_coords_file $scalar_h_coords_file $stretching_parameter $orography_layers $type_of_vertical_grid $toa
fi
fi
if [ $? -ne 0 ]
then
echo -e ${RED}Grid file creation failed.$NC
else
echo "Grid file created sucessfully."
fi
