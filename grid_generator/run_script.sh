#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

# oro_id	description
# 0			no orography
# 1			Gaussian mountain at 0 N / 0 E, H = 10 km
# 2			JW-test orography
# 3			real data interpolated to model grid
# See handbook for more information.

oro_id=1 # The orography ID. A corresponding file must exist in orography_generator/orographies.
optimize=0 # Determines wether or not the grid will be optimized.
n_iterations=2000 # The number of iterations to be used for the optimization. Only relevant, if optimize == 1.
use_scalar_h_coords_file=1 # If this is set to one, the horizontal coordinates of the grid points will be read from the file specified in the next line.
scalar_h_coords_file="grids/B5L26T30000_O0_OL17_SCVT.nc" # File used for reading horizontal coordinates of grid points, if use_scalar_h_coords_file == 1.
stretching_parameter=1.0 # stretching parameter of the vertical grid, must be >= 1, 1: no stretching
valgrind_check=0 # set this to 1, if you want to check the code with Valgrind
number_of_cpus=1 # relevant only for MPI
export OMP_NUM_THREADS=7 # relevant only for OMP
source .sh/run.sh
