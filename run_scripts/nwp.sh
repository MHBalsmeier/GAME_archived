#!/bin/bash

# basic run properties
operator=MHB
nwp_mode=0 # numerical weather prediction mode (calls data assimilation, special output, etc.)
# test_id must only be set if nwp_mode == 0
test_id=3
overwrite_run_id=1
# run_id must only be set if overwrite_run_id == 0 (otherwise it is chosen automatically)
run_id=ullrich_dry_rev
run_span=43200

# grid properties
res_id=5 #  resolution ID (number of bisections of basic icosahedral triangles)
number_of_layers=26
toa=30000 # top of atmosphere
orography_layers=17 # number of layers following the orography
# orography_id must only be set if test_id == -1. Otherwise the model will know which orography_id to use from the test_id.
orography_id=2

# "physics" configuration
cfl_margin=0.25 # The time step will be calculated as follows (delta t) = (delta t from CFL)*(1 - cfl_margin). 0.25 can be considered a safe standard value.
temperature_diff_h=0 # turn on if you want horizontal temperature diffusion
temperature_diff_v=0 # turn on if you want vetical temperature diffusion
momentum_diff=0 # turn on if you momentum diffusion
tracers_on=0 # turn on if you want advect tracers and have phase transitions
rad_on=0 # turn on if you want radiation
radiation_delta_t=3600 # every how many seconds the radiation fluxes wil be updated

# I/O
init_dir_base=input # where the program will look for the NC file containing the initialization state
output_dir_base=output # where the model will create the output directory
start_year=2000 # not relevant in NWP mode
start_month=1 # not relevant in NWP mode
start_day=1 # not relevant in NWP mode
start_hour=0 # not relevant in NWP mode
write_out_interval=900 # every how many seconds an output file will be created
write_out_mass_dry_integral=1
write_out_entropy_gas_integral=1
write_out_linearized_entropy_gas_integral=1
write_out_energy_integral=1

# parallelization
# relevant only for OMP
export OMP_NUM_THREADS=5
# relevant only for MPI
number_of_cpus=1

# data assimilation
ndvar_directory=/home/max/compiled/ndvar # The directory, where the ndvar program is stored. Only relevant in NWP mode.

# that's it, now the basic run script will be called
source core/run.sh




