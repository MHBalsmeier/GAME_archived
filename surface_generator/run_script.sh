#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# oro_id	description
# 1			Gaussian mountain at 0 N / 0 E, H = 10 km
# 2			real data interpolated to model grid
# See handbook for more information.

oro_id=0 # The orography ID.
oro_rescale_factor=0.5 # rescale factor for the orography (0.5 is sufficient to stabilize the model in real orography (oro_id=2))
res_id=5 # resolution ID
valgrind_check=0 # set this to 1, if you want to check the code with Valgrind
export OMP_NUM_THREADS=1 # relevant only for OMP
source .sh/root_script.sh
