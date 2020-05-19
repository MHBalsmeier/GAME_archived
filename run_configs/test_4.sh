#!/bin/bash
run_span=400
write_out_interval=200
grid_props_file=/home/max/compiled/game/grids/B4L6T30000_M2_O1_OL4.nc
init_state_file=/home/max/compiled/game/input/test_4_B4L6T30000_M2_O1_OL4.grb2
output_dir=/home/max/compiled/game/output/test_4
cfl_margin=0.2
dissipation=1
rad_on=1
add_comps_on=1
source run.sh
