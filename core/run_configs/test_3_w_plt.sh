#!/bin/bash
run_span=36000
write_out_interval=288
grid_props_file=/home/max/my_code/game/grid_generator/nc_files/B4L6T30000_M2_O0.nc
init_state_file=/home/max/my_code/game/test_generator/grib_files/test_3_B4L6T30000_M2_O0.grb2
output_dir=output/test_3
cfl_margin=0.2
source run.sh
disp_level=6
disp_shortname=prmsl
fig_save_path=/home/max/figs/game_output/movie_elements
core_path=/home/max/my_code/game/core
source create_movie.sh
