# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

oro_id=3
# relevant only for OMP
export OMP_NUM_THREADS=5
# relevant only for MPI
number_of_cpus=1
source .sh/run.sh
