#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

game_home_dir=/home/max/compiled/game_dev
run_id=ullrich_dry_rev
output_dir=$game_home_dir/output/$run_id
fig_save_path=/home/max/figs/game_output/$run_id
disp_shortname=2t
disp_level=2
source $game_home_dir/plotting/.sh/movie_root.sh
