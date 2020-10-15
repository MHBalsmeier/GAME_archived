#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

# basic run properties
game_home_dir=/home/max/compiled/game_dev # the root directory of your GAME instance
operator=MHB # the ID of the person / group / institution running the model
nwp_mode=0 # switches on or off the numerical weather prediction mode (calls data assimilation, special output, etc.)
test_id=3 # Must only be set if nwp_mode == 0. In this case, it specifies which test scenario to run.
overwrite_run_id=1 # this overwrites the otherwise automatically chosen run_id
run_id=ullrich_dry_irrev # run_id must only be set if overwrite_run_id == 0 (otherwise it is chosen automatically)
run_span=43200 # how long the model is supposed to run

# grid properties
res_id=5 #  resolution ID (number of bisections of basic icosahedral triangles)
number_of_layers=26
toa=30000 # top of atmosphere
orography_layers=17 # number of layers following the orography

# "physics" configuration
cfl_margin=0.25 # The time step will be calculated as follows (delta t) = (delta t from CFL)*(1 - cfl_margin). 0.25 can be considered a safe standard value.
mass_dry_diff_h=0 # turn on if you want horizontal dry mass diffusion
mass_dry_diff_v=0 # turn on if you want vetical dry mass diffusion
temperature_diff_h=1 # turn on if you want horizontal temperature diffusion
temperature_diff_v=1 # turn on if you want vetical temperature diffusion
momentum_diff=1 # turn on if you want momentum diffusion
tracers_on=0 # turn on if you want advect tracers and have phase transitions
rad_on=0 # turn on if you want radiation
radiation_delta_t=3600 # every how many seconds the radiation fluxes wil be updated

# I/O
start_year=2000 # not relevant in NWP mode
start_month=1 # not relevant in NWP mode
start_day=1 # not relevant in NWP mode
start_hour=0 # not relevant in NWP mode
write_out_interval=900 # every how many seconds an output file will be created
write_out_mass_dry_integral=1 # If set to 1, the total dry mass of the atmosphere will be written out at every time step.
write_out_entropy_gas_integral=1 # If set to 1, the total entropy of the atmosphere will be written out at every time step.
write_out_linearized_entropy_gas_integral=1 # If set to 1, the total linearized entropy (proportional to density times potential temperature) of the atmosphere will be written out at every time step.
write_out_energy_integral=1 # If set to 1, the total integrals of the energy forms of the atmosphere will be written out at every time step.
grib_output_switch=1 # If set to 1, output will be written to grib files.
netcdf_output_switch=0 # If set to 1, output will be written to netcdf files.
synop_output_mode=0 # If set to 1, additional output on synoptical pressure levels will be created. The pressure levels can be set in the file core/src/settings.c. The numer of pressure levels must be set in the file core/src/settings.h.
aviation_output_mode=0 # If set to 1, additional output on flight levels will be created. The flight levels can be set in the file core/src/settings.c. The numer of flight levels must be set in the file core/src/settings.h.

# parallelization
# relevant only for OMP
export OMP_NUM_THREADS=1
# relevant only for MPI
number_of_cpus=1

# data assimilation
ndvar_directory=/home/max/compiled/ndvar # The directory, where the ndvar program is stored. Only relevant in NWP mode.

# that's it, now the basic run script will be sourced
source $game_home_dir/core/run.sh




