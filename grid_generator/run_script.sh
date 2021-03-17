#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/AUN4GFD/game

# oro_id	description
# 0			no orography
# 1			Gaussian mountain at 0 N / 0 E, H = 10 km
# 2			JW-test orography
# 3			real data interpolated to model grid
# See handbook for more information.

oro_id=3 # The orography ID. A corresponding file must exist in orography_generator/orographies.
optimize=0 # Determines wether or not the grid will be optimized.
n_iterations=2000 # The number of iterations to be used for the optimization. Only relevant, if optimize == 1.
use_scalar_h_coords_file=1 # If this is set to one, the horizontal coordinates of the grid points will be read from the file specified in the next line.
scalar_h_coords_file="grids/B5L26T30000_O0_OL23_SCVT.nc" # File used for reading horizontal coordinates of grid points, if use_scalar_h_coords_file == 1.
stretching_parameter=1.0 # stretching parameter of the vertical grid, must be >= 1, 1: no stretching
orography_layers=25 # number of layers following orography (only relevant if type_of_vertical_grid == 0)
toa=41152 # height of the top of the atmosphere
type_of_vertical_grid=0 # 0: terrain following coordinates, 1: block-like orography
valgrind_check=0 # set this to 1, if you want to check the code with Valgrind
number_of_cpus=1 # relevant only for MPI
export OMP_NUM_THREADS=6 # relevant only for OMP
source .sh/run.sh
