#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This is a run script for operational runs (NWP runs).

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
game_home_dir=${BASH_ARGV[4]}
ideal_input_id=-1 # specifies which test scenario to run (-1 corresponds to an NWP run)
run_span=${BASH_ARGV[5]} # how long the model is supposed to run; for small Earth experiments this will be rescaled proportional to the radius
run_id=${BASH_ARGV[6]} # how long the model is supposed to run
valgrind_check=0 # set this to 1, if you want to check the code with Valgrind
start_year=${BASH_ARGV[3]} # defines the start time of the model run
start_month=${BASH_ARGV[2]} # defines the start time of the model run
start_day=${BASH_ARGV[1]} # defines the start time of the model run
start_hour=${BASH_ARGV[0]} # defines the start time of the model run
orography_id=${BASH_ARGV[7]} # ID of the orography field. Based on this the grid file will be chosen.

# diffusion settings
momentum_diff_h=1 # turn on if you want horizontal momentum diffusion
momentum_diff_v=1 # turn on if you want vertical momentum diffusion
temperature_diff_h=1 # turn on if you want horizontal temperature diffusion
temperature_diff_v=1 # turn on if you want vetical temperature diffusion
tracer_diff_h=1 # turn on if you want horizontal tracer concentration diffusion
tracer_diff_v=1 # turn on if you want vertical tracer concentration diffusion

# "physics" configuration
assume_lte=1 # set this to one if you do not want to assign individual temperatures to tracers
rad_on=1 # set to 0 if you want no radiation, 1 for real radiation and 2 for Held-Suarez forcing
no_rad_moisture_layers=12 # number of layers in which radiation-moisture interaction is neglected
soil_on=1 # switch for the soil component of the model
explicit_boundary_layer=0 # switch for an additional simplified horizontal friction in the boundary layer

# I/O
write_out_interval=10800 # every how many seconds an output file will be created; for small Earth experiments this will be rescaled proportional to the radius
write_out_integrals=0 # If set to 1, fundamental integrals of the atmosphere will be written out at every time step.
model_level_output_switch=0 # If set to 1, variables will be written out on model levels.
pressure_level_output_switch=1 # If set to 1, additional output on pressure levels will be created. The pressure levels can be set in the file src/io/write_output.c.
surface_output_switch=1 # If set to 1, surface variables will be diagnozed and writing to separate files.
grib_output_switch=1 # If set to 1, output will be written to grib files.
netcdf_output_switch=0 # If set to 1, output will be written to netcdf files.
time_to_next_analysis=${BASH_ARGV[8]} # the time between this model run and the next analysis, only relevant in NWP mode for data assimilation

# parallelization
export OMP_NUM_THREADS=${BASH_ARGV[9]} # relevant for OMP

# that's it, now the basic run script will be sourced
source $game_home_dir/run_scripts/.sh/root_script.sh




