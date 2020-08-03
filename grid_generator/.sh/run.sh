# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game


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
echo "Horizontal coordinates of the generating points (the scalar points in terms of the model) will be read from file $scalar_h_coords_file."
fi
echo ""
echo "********** Calling the GAME grid generator **********"
echo ""
mpirun -np $number_of_cpus ./grid_generator $oro_id $optimize $n_iterations $use_scalar_h_coords_file $scalar_h_coords_file
if [ $? -ne 0 ]
then
echo -e ${RED}Grid file creation failed.$NC
else
echo "Grid file created sucessfully."
fi
