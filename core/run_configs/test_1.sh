#!/bin/bash
run_span=7200
write_out_interval=100
geo_prop_file=../grid_generator/nc_files/res_4_oro_0_geo_prop.nc
init_state_file=../test_generator/grib_files/test_1_res_4_oro_0.grb2
output_folder=output/test_1
source run.sh
