#!/bin/bash

# This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

# oroid		description
# 0			no orography
# 1			Gaussian mountain at 0 N / 0 E, H = 10 km
# 2			JW-test orography
# 3			real data interpolated to model grid
# See handbook for more information.

oro_id=2
optimize=0
n_iterations=2000
use_scalar_h_coords_file=1
scalar_h_coords_file="grids/B6L26T30000_O0_OL17_SCVT.nc"
# relevant only for MPI
number_of_cpus=1
# relevant only for OMP
export OMP_NUM_THREADS=5
source .sh/run.sh
