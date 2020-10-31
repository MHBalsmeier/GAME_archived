#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

game_home_dir=${BASH_ARGV[0]} # the home directory of GAME
run_id=${BASH_ARGV[1]} # the run id which you want to plot
run_span=${BASH_ARGV[2]} # the length of the run
output_dir=$game_home_dir/output/$run_id # the directory where the grib files are stored
fig_save_path=~/figs/game_output/$run_id # the path to which the figures will be saved
grid_props_file=$game_home_dir/grids/B5L26T30000_O0_OL17_SCVT.nc # the file where the grid properties are stored
disp_shortname_list=(2t prmsl t r gh 2t prmsl t r gh) # short names according to grib as an array 
disp_level_list=(2 0 850 850 500 2 0 850 850 500) # levels according to grib as an array
on_pressure_level_list=(0 0 1 1 1 0 0 1 1 1) # set this to 1 for each plot individually if the variable resides on pressure levels
plot_intervals_list=(10800 21600 10800 10800 21600 10800 21600 10800 10800 21600) # every how many seconds you want to plot each variable
uniform_colormap_list=(0 0 0 0 0 0 0 0 0 0) # set this to 1 for each plot individually if you want to enforce a uniform colormap for all the time steps
scope_list=(CONUS CONUS CONUS CONUS CONUS CEU CEU CEU CEU CEU) # the areas of the plots
projections_list=(Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic) # the projections of the plots
source $game_home_dir/plotting/.sh/maps_root.sh # this is the script from which the python plot scripts are called
