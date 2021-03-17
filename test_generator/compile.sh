#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/AUN4GFD/game

dcmip_home_dir=../../DCMIP2016 # set this to the directory of the DCMIP2016 repository

echo "Starting to compile test generator."
mpicc src/test_generator.c ../core/src/spatial_operators/inner_product.c ../core/src/spatial_operators/gradient_operators.c ../core/src/spatial_operators/scalar_times_vector.c ../core/src/spatial_operators/vorticities.c ../core/src/spatial_operators/vorticity_flux.c ../core/src/diagnostics/averaging.c ../core/src/diagnostics/thermodynamic_diagnostics.c ../core/src/settings.c ../core/src/io/set_grid_props_and_dt.c $dcmip_home_dir/interface/baroclinic_wave_test.f90 ../shared/various.c -leccodes -lnetcdf -lm -lgeos95 -fopenmp -latmostracers -Wl,-rpath=/lib -Wall -o test_generator
if [ $? -ne 0 ]
then
echo -e ${RED}Test generator compilation failed.$NC
else
echo "Test generator compiled sucessfully."
fi
