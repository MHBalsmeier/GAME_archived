#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This is a run script for the Held-Suarez testcase.

# test_id	description
# 0			standard atmosphere
# 1			dry Ullrich test
# 2			moist Ullrich test
# -1 	  	NWP run
# oro_id	description
# 0			no orography
# 1			real data interpolated to the model grid
# See handbook for more information.

# basic run properties
game_home_dir=/home/max/code/GAME
ideal_input_id=0 # specifies which test scenario to run
run_id=held_suar # run_id must only be set if ideal_input_id != -1 (otherwise it is chosen automatically)
run_span=$((1250*24*3600)) # how long the model is supposed to run; for small Earth experiments this will be rescaled proportional to the radius
valgrind_check=0 # set this to 1 if you want to check the code with Valgrind
start_year=2000 # defines the start time of the model run
start_month=1 # defines the start time of the model run
start_day=1 # defines the start time of the model run
start_hour=0 # defines the start time of the model run

# grid properties
orography_id=0 # ID of the orography field. Based on this the grid file will be chosen.

# dynamics settings
dt_parameter=1.5 # The sound time step will be calculated as follows: delta_t = dt_parameter*delta_x / km. 1.5 can be considered a standard value.
momentum_diff_h=1 # turn on if you want horizontal momentum diffusion
momentum_diff_v=0 # turn on if you want vertical momentum diffusion
diff_h_smag_div=0.01 # horizontal diffusion Smagorinsky factor acting on divergent movements
diff_h_smag_rot=0.01 # horizontal diffusion Smagorinsky factor acting on vortical movements
temperature_diff_h=0 # turn on if you want horizontal temperature diffusion
temperature_diff_v=0 # turn on if you want vetical temperature diffusion
tracer_diff_h=0 # turn on if you want horizontal tracer concentration diffusion
tracer_diff_v=0 # turn on if you want vertical tracer concentration diffusion
damping_start_height_over_toa=0.53 # swamp layer boundary in relation to the TOA
damping_coeff_max=0.25 # maximum swamp layer damping coefficient
explicit_boundary_layer=1 # switch for an additional simplified horizontal friction in the boundary layer

# "physics" configuration
rad_on=2 # set to 0 if you want no radiation, 1 for real radiation and 2 for Held-Suarez forcing
radiation_delta_t=3600 # every how many seconds the radiation fluxes wil be updated; for small Earth experiments this will be rescaled proportional to the radius
assume_lte=1 # set this to one if you do not want to assign individual temperatures to tracers

# I/O
write_out_interval=86400 # every how many seconds an output file will be created; for small Earth experiments this will be rescaled proportional to the radius
write_out_integrals=1 # If set to 1, fundamental integrals of the atmosphere will be written out at every time step.
model_level_output_switch=1 # If set to 1, variables will be written out on model levels.
pressure_level_output_switch=1 # If set to 1, additional output on pressure_leveltical pressure levels will be created. The pressure levels can be set in the file core/src/settings.c. The number of pressure levels must be set in the file core/src/settings.h.
surface_output_switch=1 # If set to 1, surface variables will be diagnozed and writing to separate files.
grib_output_switch=1 # If set to 1, output will be written to grib files on a lat-lon grid.
netcdf_output_switch=0 # If set to 1, output will be written to netcdf files on the hexagonal (and pentagonal) cell centers.
delta_t_between_analyses=-1 # the time difference between two analyses, only relevant in NWP mode

# parallelization
export OMP_NUM_THREADS=2 # relevant for OMP

# that's it, now the basic run script will be sourced
source $game_home_dir/run_scripts/.sh/root_script.sh




