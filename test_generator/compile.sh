#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

dcmip_home_dir=../../DCMIP2016 # set this to the directory of the DCMIP2016 repository

echo "Starting to compile test generator."
gcc -O2 src/test_generator.c ../src/spatial_operators/inner_product.c ../src/spatial_operators/gradient_operators.c ../src/spatial_operators/multiplications.c ../src/spatial_operators/vorticities.c ../src/spatial_operators/vorticity_flux.c ../src/spatial_operators/averaging.c ../src/thermodynamics.c ../src/io/set_grid_props_and_dt.c ../grid_generator/src/vertical_grid.c $dcmip_home_dir/interface/baroclinic_wave_test.f90 -leccodes -lnetcdf -lm -lgeos95 -fopenmp -latmostracers -Wall -o test_generator
if [ $? -ne 0 ]
then
echo -e ${RED}Test generator compilation failed.$NC
else
echo "Test generator compiled successfully."
fi
