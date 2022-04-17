#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

run_dir=$game_home_dir/output/$run_id
if [ -d $run_dir ]
then
rm -r $run_dir
fi
mkdir $run_dir
cd $run_dir

if [ ! -f $game_home_dir/build/game ]
then
echo "Executable game missing. Compile first. Aborting run."
cd - > /dev/null
exit 1
fi

if [ -f game ]
then
rm game
fi

cp $game_home_dir/build/game .

if [ $valgrind_check -eq 0 ]
then
./game $run_span $write_out_interval $momentum_diff_h $momentum_diff_v $rad_on $soil_on $no_rad_moisture_layers $write_out_integrals $temperature_diff_h $start_year $start_month $start_day $start_hour $temperature_diff_v $run_id $orography_id $ideal_input_id $grib_output_switch $netcdf_output_switch $pressure_level_output_switch $model_level_output_switch $surface_output_switch $assume_lte $delta_t_between_analyses $damping_start_height_over_toa $damping_coeff_max $explicit_boundary_layer $tracer_diff_h $tracer_diff_v
fi
if [ $valgrind_check -eq 1 ]
then
valgrind ./game $run_span $write_out_interval $momentum_diff_h $momentum_diff_v $rad_on $soil_on $no_rad_moisture_layers $write_out_integrals $temperature_diff_h $start_year $start_month $start_day $start_hour $temperature_diff_v $run_id $orography_id $ideal_input_id $grib_output_switch $netcdf_output_switch $pressure_level_output_switch $model_level_output_switch $surface_output_switch $assume_lte $delta_t_between_analyses $damping_start_height_over_toa $damping_coeff_max $explicit_boundary_layer $tracer_diff_h $tracer_diff_v
fi
cd - > /dev/null
