#!/bin/bash

run_id=ullrich_dry_rev
run_span=43200
output_dir=/home/max/compiled/game_dev/output/$run_id
grid_props_file=
disp_shortname_list=(2t prmsl 2t prmsl)
disp_level_list=(2 0 2 0)
plot_intervals_list=(10800 21600 10800 21600)
uniform_colormap_list=(0 0 0 0)
scope_list=(CONUS CONUS CEU CEU)
projections_list=(Gnomonic Gnomonic Gnomonic Gnomonic)
fig_save_path=~/figs/game_output/$run_id
grid_props_file=/home/max/compiled/game_dev/grids/B5L26T30000_O2_OL17_SCVT.nc
source plotting/.sh/maps_root.sh
