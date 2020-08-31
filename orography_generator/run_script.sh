#!/bin/bash

# This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

# oroid		description
# 1			Gaussian mountain at 0 N / 0 E, H = 10 km
# 2			JW-test orography
# 3			real data interpolated to model grid
# See handbook for more information.

oro_id=3
# relevant only for OMP
export OMP_NUM_THREADS=5
# relevant only for MPI
number_of_cpus=1
source .sh/run.sh
