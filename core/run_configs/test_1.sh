#!/bin/bash
run_span=1800
write_out_interval=288
grid_props_file=/home/max/my_code/game/grid_generator/nc_files/B4L6T30000_M2_O0.nc
init_state_file=/home/max/my_code/game/test_generator/grib_files/test_1_B4L6T30000_M2_O0.grb2
output_dir=output/test_1
cfl_margin=0.2
source run.sh
