#!/bin/bash
run_span=1800
write_out_interval=100
grid_props_file=/home/max/my_code/game/grid_generator/nc_files/B4L6T20000_O0.nc
init_state_file=/home/max/my_code/game/test_generator/grib_files/test_3_B4L6T20000_O0.grb2
output_dir=output/test_3
source run.sh
disp_level=2
disp_shortname=u
fig_save_path=/home/max/figs/game_output/movie_elements
core_path=/home/max/my_code/game/core
python /home/max/my_code/game_plotting/create_movie_elements.py $run_span $write_out_interval $disp_level $disp_shortname $grid_props_file $fig_save_path $core_path $output_dir
if [ -f movie.mp4 ]
then
rm movie.mp4
fi
/home/max/my_code/game_plotting/create_movie.sh
