#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# oro_id	description
# 1			Gaussian mountain at 0 N / 0 E, H = 10 km
# 2			real data interpolated to model grid
# See handbook for more information.

oro_id=2 # The orography ID.
no_of_cells_for_gliding_avrg=9 # number of cells of the model grid used for averaging the input dataset to the scalar data points (increase for smoothing)
res_id=5 # resolution ID
valgrind_check=0 # set this to 1, if you want to check the code with Valgrind
export OMP_NUM_THREADS=1 # relevant only for OMP
source .sh/root_script.sh
