#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

# oroid		description
# 1			Gaussian mountain at 0 N / 0 E, H = 10 km
# 2			JW-test orography
# 3			real data interpolated to model grid
# See handbook for more information.

oro_id=2 # The orography ID.
valgrind_check=0 # set this to 1, if you want to check the code with Valgrind
export OMP_NUM_THREADS=5 # relevant only for OMP
number_of_cpus=1 # relevant only for MPI
source .sh/run.sh
