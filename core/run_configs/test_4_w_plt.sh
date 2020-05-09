#!/bin/bash
run_span=120000
write_out_interval=12000
grid_props_file=/home/max/compiled/game/grids/B4L6T30000_M2_O0.nc
init_state_file=/home/max/compiled/game/input/test_4_B4L6T30000_M2_O0.grb2
output_dir=/home/max/compiled/game/output/test_4
cfl_margin=0.2
dissipation=1
rad_on=1
add_comps_on=1
source run
disp_level=2
disp_shortname=2t
fig_save_path=/home/max/figs/game_output/movie_elements
core_path=/home/max/my_code/game/core
source create_movie.sh
