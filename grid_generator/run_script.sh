# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

oro_id=2
optimize=1
n_iterations=0
use_scalar_h_coords_file=1
scalar_h_coords_file="nc_files/B5L26T30000_O0_OL17_SCVT.nc"
# relevant only for MPI
number_of_cpus=1
# relevant only for OMP
export OMP_NUM_THREADS=5
source .sh/run.sh
