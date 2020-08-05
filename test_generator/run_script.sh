# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

# test_ids:
# 0:	standard atmosphere without orography
# 1:	standard atmosphere with Gaussian mountain
# 2:	JW test, dry, balanced
# 3:	JW test, dry, perturbed
# 4:	JW test, moist, balanced
# 5:	JW test, moist, perturbed

test_id=5
# relevant only for MPI
number_of_cpus=1
# relevant only for OMP
export OMP_NUM_THREADS=5
source .sh/run.sh
