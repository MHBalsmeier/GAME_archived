#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

game_home_dir=/home/max/compiled/game # the home directory of GAME
run_id=ideal # the run id which you want to plot
run_span=$((8*86400)) # the length of the run
output_dir=$game_home_dir/output/$run_id # the directory where the grib files are stored
fig_save_path=/home/max/figs/game_output # the path to which the figures will be saved
init_year=2000 # year of the start of the model run
init_month=1 # month of the start of the model run
init_day=1 # day of the start of the model run
init_hr=0 # hour of the start of the model run
start_time_since_init=0 # when to begin plotting reative to the model initialization
disp_shortname_list=(
sp) # short names according to grib as an array 
disp_level_list=(
0) # levels according to grib as an array
on_pressure_level_list=(
0) # set this to 1 for each plot individually if the variable resides on pressure levels
plot_intervals_list=(
86400) # every how many seconds you want to plot each variable
uniform_colormap_list=(
1) # set this to 1 for each plot individually if you want to enforce a uniform colormap for all the time steps
scope_list=(
WORLD) # the areas of the plots
projections_list=(
EckertIII) # the projections of the plots
synoptical_time_mode=(
1) # this forces the time description to be of the form "init: ..., valid: ... (+ ....)"
source $game_home_dir/plotting/.sh/maps_root.sh # this is the script from which the python plot scripts are called
