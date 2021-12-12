#!/bin/bash

game_home_dir=${BASH_ARGV[2]} # the home directory of GAME
run_id=${BASH_ARGV[1]} # the run id which you want to plot
run_span=${BASH_ARGV[0]} # the length of the run
output_dir=$game_home_dir/output/$run_id # the directory where the grib files are stored
fig_save_path=${BASH_ARGV[3]} # the path to which the figures will be saved
init_year=${BASH_ARGV[4]} # year of the start of the model run
init_month=${BASH_ARGV[5]} # month of the start of the model run
init_day=${BASH_ARGV[6]} # day of the start of the model run
init_hr=${BASH_ARGV[7]} # hour of the start of the model run
plot_interval=${BASH_ARGV[8]} # the interval between plots in seconds
start_time_since_init=${BASH_ARGV[9]} # when to begin plotting reative to the model initialization
omp_num_threads=${BASH_ARGV[10]} # relevant only for OMP
disp_shortname_list=(
2t gust rprate sprate tcc
) # short names according to grib as an array 
disp_level_list=(
2 10 0 0 0
) # levels according to grib as an array
on_pressure_level_list=(
0 0 0 0 0
) # set this to 1 for each plot individually if the variable resides on pressure levels
plot_intervals_list=(
$plot_interval $plot_interval $plot_interval $plot_interval $plot_interval
) # every how many seconds you want to plot each variable
uniform_colormap_list=(
0 0 0 0 0
) # set this to 1 for each plot individually if you want to enforce a uniform colormap for all the time steps
scope_list=(
CEU CEU CEU CEU CEU
) # the areas of the plots
projections_list=(
Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic
) # the projections of the plots
synoptical_time_mode=(
1 1 1 1 1
) # this forces the time description to be of the form "init: ..., valid: ... (+ ....)"
source $game_home_dir/plotting/.sh/maps_root.sh # this is the script from which the python plot scripts are called




