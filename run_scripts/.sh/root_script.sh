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
./game $run_span $write_out_interval $dt_parameter $momentum_diff_h $momentum_diff_v $rad_on $write_out_mass_integrals $write_out_rhotheta_integral $write_out_energy_integral $temperature_diff_h $radiation_delta_t $start_year $start_month $start_day $start_hour $temperature_diff_v $run_id $toa $orography_id $ideal_input_id $grib_output_switch $netcdf_output_switch $pressure_level_output_switch $model_level_output_switch $surface_output_switch $orography_layers $type_of_vertical_grid $assume_lte $slow_fast_ratio $delta_t_between_analyses $diff_h_smag_div $diff_h_smag_rot $damping_start_height_over_toa $damping_coeff_max $explicit_boundary_layer $tracer_diff_h $tracer_diff_v $impl_thermo_weight $cloud_droplets_velocity $precipitation_droplets_velocity $mixing_length
fi
if [ $valgrind_check -eq 1 ]
then
valgrind ./game $run_span $write_out_interval $dt_parameter $momentum_diff_h $momentum_diff_v $rad_on $write_out_mass_integrals $write_out_rhotheta_integral $write_out_energy_integral $temperature_diff_h $radiation_delta_t $start_year $start_month $start_day $start_hour $temperature_diff_v $run_id $toa $orography_id $ideal_input_id $grib_output_switch $netcdf_output_switch $pressure_level_output_switch $model_level_output_switch $surface_output_switch $orography_layers $type_of_vertical_grid $assume_lte $slow_fast_ratio $delta_t_between_analyses $diff_h_smag_div $diff_h_smag_rot $damping_start_height_over_toa $damping_coeff_max $explicit_boundary_layer $tracer_diff_h $tracer_diff_v $impl_thermo_weight $cloud_droplets_velocity $precipitation_droplets_velocity $mixing_length
fi
cd - > /dev/null
