#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/AUN4GFD/game

game_home_dir=/home/max/compiled/game_dev # the home directory of GAME
run_id=ullrich_dry_flat # the run id which you want to plot
run_span=0 # the length of the run
output_dir=$game_home_dir/output/$run_id # the directory where the grib files are stored
fig_save_path=/home/max/figs/game_output # the path to which the figures will be saved
init_year=2000 # year of the start of the model run
init_month=1 # month of the start of the model run
init_day=1 # day of the start of the model run
init_hr=0 # hour of the start of the model run
grid_props_file=$game_home_dir/grids/B5L26T30000_O0_OL17_SCVT.nc # the file where the grid properties are stored
disp_shortname_list=(
surface_wind 2t rprate sprate gust prmsl gh gh cape tcc r) # short names according to grib as an array 
disp_level_list=(
10 2 0 0 10 0 500 200 0 0 850) # levels according to grib as an array
on_pressure_level_list=(
0 0 0 0 0 0 1 1 0 0 1) # set this to 1 for each plot individually if the variable resides on pressure levels
plot_intervals_list=(
21600 21600 21600 21600 21600 21600 21600 21600 21600 21600 21600) # every how many seconds you want to plot each variable
uniform_colormap_list=(
0 0 0 0 0 0 0 0 0 0 0) # set this to 1 for each plot individually if you want to enforce a uniform colormap for all the time steps
scope_list=(
CHINA CHINA CHINA CHINA CHINA CHINA CHINA CHINA CHINA CHINA CHINA) # the areas of the plots
projections_list=(
Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic) # the projections of the plots
synoptical_time_mode=(
1 1 1 1 1 1 1 1 1 1 1) # this forces the time description to be of the form "init: ..., valid: ... (+ ....)"
source $game_home_dir/plotting/.sh/maps_root.sh # this is the script from which the python plot scripts are called
