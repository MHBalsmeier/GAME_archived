#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/AUN4GFD/game

# basic run properties
game_home_dir=/home/max/compiled/game_dev # the root directory of your GAME installation
operator=MHB # the ID of the person / group / institution running the model
ideal_input_id=0 # specifies which test scenario to run
run_id=standard # run_id must only be set if ideal_input_id != -1 (otherwise it is chosen automatically)
run_span=$((10*24*3600)) # how long the model is supposed to run
valgrind_check=0 # set this to 1, if you want to check the code with Valgrind
start_year=2000 # defines the start time of the model run
start_month=1 # defines the start time of the model run
start_day=1 # defines the start time of the model run
start_hour=0 # defines the start time of the model run

# grid properties
toa=41152 # top of atmosphere
type_of_vertical_grid=0 # 0: terrain following coordinates, 1: block-like orography
orography_id=0 # ID of the orography field, for ideal_input_id > -1 orography_id will be set automatically
orography_layers=23 # number of layers following orography (only relevent if type_of_vertical_grid == 0)

# dynamics settings
cfl_margin=0.55 # The sound time step will be calculated as follows (delta t) = (delta t from horizontal CFL)*(1 - cfl_margin). 0.55 can be considered a safe standard value.
adv_sound_ratio=1 # the ratio of the advective to the sound time step
momentum_diff_h=1 # turn on if you want horizontal momentum diffusion
dissipative_heating=1 # turn on if you want a dissipative heating rate
diff_h_smag_fac=0.13 # horizontal diffusion Smagorinsky factor
shear_bg=1.5e-5 # assumed background (minimum) shear
mass_advection_order=3 # set this to 2 or 3 (convergence order of the mass advection)
entropy_advection_order=3 # set this to 2 or 3 (convergence order of the entropy advection)

# "physics" configuration
mass_dry_diff_h=0 # turn on if you want horizontal dry mass diffusion
mass_dry_diff_v=0 # turn on if you want vetical dry mass diffusion
temperature_diff_h=0 # turn on if you want horizontal temperature diffusion
temperature_diff_v=0 # turn on if you want vetical temperature diffusion
momentum_diff_v=0 # turn on if you want vertical momentum diffusion
damping_start_height_over_toa=0.53 # Swamp layer boundary in relation to the TOA.
damping_coeff_max=0.25 # maximum swamp layer damping coefficient
rad_on=0 # turn on if you want radiation
radiation_delta_t=10800 # every how many seconds the radiation fluxes wil be updated
assume_lte=1 # set this to one if you do not want to assign individual temperatures to tracers

# I/O
write_out_interval=86400 # every how many seconds an output file will be created
write_out_mass_dry_integral=1 # If set to 1, the total dry mass of the atmosphere will be written out at every time step.
write_out_entropy_gas_integral=1 # If set to 1, the total entropy of the atmosphere will be written out at every time step.
write_out_linearized_entropy_gas_integral=1 # If set to 1, the total linearized entropy (proportional to density times potential temperature) of the atmosphere will be written out at every time step.
write_out_energy_integral=1 # If set to 1, the total integrals of the energy forms of the atmosphere will be written out at every time step.
model_level_output_switch=0 # If set to 1, variables will be written out on model levels.
pressure_level_output_switch=1 # If set to 1, additional output on pressure_leveltical pressure levels will be created. The pressure levels can be set in the file core/src/settings.c. The numer of pressure levels must be set in the file core/src/settings.h.
surface_output_switch=1 # If set to 1, surface variables will be diagnozed and writing to separate files.
grib_output_switch=1 # If set to 1, output will be written to grib files on a lat-lon grid.
netcdf_output_switch=0 # If set to 1, output will be written to netcdf files on the hexagonal (and pentagonal) cell centers.
delta_t_between_analyses=-1 # the time difference between two analyses, only relevant in NWP mode

# parallelization
export OMP_NUM_THREADS=6 # relevant only for OMP
number_of_cpus=1 # relevant only for MPI

# that's it, now the basic run script will be sourced
source $game_home_dir/core/run.sh




