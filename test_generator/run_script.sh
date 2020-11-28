#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

# test_ids:
# 0:	standard atmosphere without orography
# 1:	standard atmosphere with Gaussian mountain
# 2:	JW test, dry, balanced
# 3:	JW test, dry, perturbed
# 4:	JW test, moist, balanced
# 5:	JW test, moist, perturbed
# 6:	JW test, dry, balanced, with oro_id = 3
# 7:	JW test, moist, balanced, with oro_id = 3
# 8:	Ullrich test, dry, perturbed
# 9:	Ullrich test, moist, perturbed
# 10:	Ullrich test, dry, perturbed, with oro_id = 3
# 11:	Ullrich test, moist, perturbed, with oro_id = 3
# 12:	standard atmosphere with oro_id = 3

test_id=1
orography_layers=17 # number of layers following orography
valgrind_check=0 # set this to 1, if you want to check the code with Valgrind
number_of_cpus=1 # relevant only for MPI
export OMP_NUM_THREADS=5 # relevant only for OMP
source .sh/run.sh
