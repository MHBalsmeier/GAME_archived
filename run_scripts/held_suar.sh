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
start_year=2000 # defines the start time of the model run
start_month=1 # defines the start time of the model run
start_day=1 # defines the start time of the model run
start_hour=0 # defines the start time of the model run
orography_id=0 # ID of the orography field. Based on this the grid file will be chosen.

# diffusion settings
momentum_diff_h=1 # turn on if you want horizontal momentum diffusion
momentum_diff_v=0 # turn on if you want vertical momentum diffusion
temperature_diff_h=0 # turn on if you want horizontal temperature diffusion
temperature_diff_v=0 # turn on if you want vetical temperature diffusion
tracer_diff_h=0 # turn on if you want horizontal tracer concentration diffusion
tracer_diff_v=0 # turn on if you want vertical tracer concentration diffusion

# "physics" configuration
assume_lte=1 # set this to one if you do not want to assign individual temperatures to tracers
rad_on=2 # set to 0 if you want no radiation, 1 for real radiation and 2 for Held-Suarez forcing
no_rad_moisture_layers=12 # number of layers in which radiation-moisture interaction is neglected
prog_soil_temp=0 # switch for prognostic soil temperature
sfc_phase_trans=0 # switch for phase transitions at the surface
sfc_sensible_heat_flux=0 # switch for sensible heat flux at the surface
pbl_scheme=2 # planetary boundary scheme: 0: off, 1: NWP, 2: Held-Suarez

# I/O
write_out_interval=86400 # every how many seconds an output file will be created; for small Earth experiments this will be rescaled proportional to the radius
write_out_integrals=1 # If set to 1, fundamental integrals of the atmosphere will be written out at every time step.
model_level_output_switch=1 # If set to 1, variables will be written out on model levels.
pressure_level_output_switch=1 # If set to 1, additional output on pressure levels will be created. The pressure levels can be set in the file src/io/write_output.c.
surface_output_switch=1 # If set to 1, surface variables will be diagnozed and writing to separate files.
grib_output_switch=1 # If set to 1, output will be written to grib files on a lat-lon grid.
netcdf_output_switch=0 # If set to 1, output will be written to netcdf files on the hexagonal (and pentagonal) cell centers.
time_to_next_analysis=-1 # the time between this model run and the next analysis, only relevant in NWP runs for data assimilation

# parallelization
export OMP_NUM_THREADS=2 # relevant for OMP

# that's it, now the basic run script will be sourced
source $game_home_dir/run_scripts/.sh/root_script.sh




