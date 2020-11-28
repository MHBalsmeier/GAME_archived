#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

# basic run properties
game_home_dir=/home/max/compiled/game_dev # the root directory of your GAME instance
operator=MHB # the ID of the person / group / institution running the model
ideal_input_id=8 # specifies which test scenario to run
run_id=ullrich_dry_flat # run_id must only be set if ideal_input_id != -1 (otherwise it is chosen automatically)
run_span=777600 # how long the model is supposed to run
valgrind_check=0 # set this to 1, if you want to check the code with Valgrind

# grid properties
toa=30000 # top of atmosphere
orography_layers=23 # number of layers following orography

# "physics" configuration
cfl_margin=0.25 # The time step will be calculated as follows (delta t) = (delta t from CFL)*(1 - cfl_margin). 0.25 can be considered a safe standard value.
mass_dry_diff_h=0 # turn on if you want horizontal dry mass diffusion
mass_dry_diff_v=0 # turn on if you want vetical dry mass diffusion
temperature_diff_h=0 # turn on if you want horizontal temperature diffusion
temperature_diff_v=0 # turn on if you want vetical temperature diffusion
momentum_diff=0 # turn on if you want momentum diffusion
rad_on=0 # turn on if you want radiation
radiation_delta_t=3600 # every how many seconds the radiation fluxes wil be updated

# I/O
start_year=2000 # defines the start time of the model run
start_month=1 # defines the start time of the model run
start_day=1 # defines the start time of the model run
start_hour=0 # defines the start time of the model run
write_out_interval=900 # every how many seconds an output file will be created
write_out_mass_dry_integral=1 # If set to 1, the total dry mass of the atmosphere will be written out at every time step.
write_out_entropy_gas_integral=1 # If set to 1, the total entropy of the atmosphere will be written out at every time step.
write_out_linearized_entropy_gas_integral=1 # If set to 1, the total linearized entropy (proportional to density times potential temperature) of the atmosphere will be written out at every time step.
write_out_energy_integral=1 # If set to 1, the total integrals of the energy forms of the atmosphere will be written out at every time step.
model_level_output_switch=0 # If set to 1, variables will be written out on model levels.
pressure_level_output_switch=1 # If set to 1, additional output on pressure_leveltical pressure levels will be created. The pressure levels can be set in the file core/src/settings.c. The numer of pressure levels must be set in the file core/src/settings.h.
surface_output_switch=1 # If set to 1, surface variables will be diagnozed and writing to separate files.
flight_level_output_switch=0 # If set to 1, additional output on flight levels will be created. The flight levels can be set in the file core/src/settings.c. The numer of flight levels must be set in the file core/src/settings.h.
grib_output_switch=1 # If set to 1, output will be written to grib files.
netcdf_output_switch=0 # If set to 1, output will be written to netcdf files.

# parallelization
export OMP_NUM_THREADS=1 # relevant only for OMP
number_of_cpus=1 # relevant only for MPI

# that's it, now the basic run script will be sourced
source $game_home_dir/core/run.sh




