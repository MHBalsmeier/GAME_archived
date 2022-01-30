#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# oro_id	description
# 0			no orography
# 1			real data interpolated to the model grid
# See handbook for more information.

res_id=5 # the resolution ID (also needs to be set in src/game_types.h)
oro_id=0 # The orography ID. A corresponding file must exist in surface_generator/surface_files.
n_iterations=2000 # The number of iterations to be used for the optimization.
use_scalar_h_coords_file=0 # If this is set to one, the horizontal coordinates of the grid points will be read from the file specified in the next line.
scalar_h_coords_file="grids/RES${res_id}_L26_ORO0.nc" # File used for reading horizontal coordinates of grid points, if use_scalar_h_coords_file == 1.
stretching_parameter=1.3 # stretching parameter of the vertical grid, must be >= 1, 1: no stretching
orography_layers=23 # number of layers following orography (only relevant if type_of_vertical_grid == 0)
toa=41152 # height of the top of the atmosphere
type_of_vertical_grid=0 # 0: terrain following coordinates, 1: block-like orography
valgrind_check=0 # set this to 1 if you want to check the code with Valgrind
export OMP_NUM_THREADS=6 # relevant only for OMP
source .sh/run.sh
