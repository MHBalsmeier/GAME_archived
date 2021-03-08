#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/AUN4GFD/game

time_string=$(date --utc +%Y%m%d%H%M%S)

# this is the ideal input case
if [ $ideal_input_id -gt "-1" ]
then
orography_id=-1 # orography_id will be set automatically, depending on the test ID
# this is the NWP case
else
orography_id=${BASH_ARGV[10]} # the orography ID can be set by the user in this case
fi

# settng the delta between analyses if it is not defined in the run script
if [ -z $delta_t_between_analyses ]
then
delta_t_between_analyses=-1
fi

output_dir=$game_home_dir/output/$run_id
if [ -d $output_dir ]
then
rm -r $output_dir
fi
mkdir $output_dir
cd $game_home_dir
if [ $valgrind_check -eq 0 ]
then
mpirun -np $number_of_cpus $game_home_dir/core/game $run_span $write_out_interval $cfl_margin $momentum_diff_h $momentum_diff_v $rad_on $operator $write_out_mass_dry_integral $write_out_entropy_gas_integral $write_out_energy_integral $temperature_diff_h $radiation_delta_t $start_year $start_month $start_day $start_hour $temperature_diff_v $run_id $write_out_linearized_entropy_gas_integral $toa $orography_id $ideal_input_id $mass_dry_diff_h $mass_dry_diff_v $grib_output_switch $netcdf_output_switch $pressure_level_output_switch $model_level_output_switch $surface_output_switch $orography_layers $type_of_vertical_grid $assume_lte $adv_sound_ratio $delta_t_between_analyses $div_damp_4th_order_switch $dissipative_heating $diff_h_smag_fac $shear_bg $damping_start_height_over_toa $damping_coeff_max
fi
if [ $valgrind_check -eq 1 ]
then
valgrind $game_home_dir/core/game $run_span $write_out_interval $cfl_margin $momentum_diff_h $momentum_diff_v $rad_on $operator $write_out_mass_dry_integral $write_out_entropy_gas_integral $write_out_energy_integral $temperature_diff_h $radiation_delta_t $start_year $start_month $start_day $start_hour $temperature_diff_v $run_id $write_out_linearized_entropy_gas_integral $toa $orography_id $ideal_input_id $mass_dry_diff_h $mass_dry_diff_v $grib_output_switch $netcdf_output_switch $pressure_level_output_switch $model_level_output_switch $surface_output_switch $orography_layers $type_of_vertical_grid $assume_lte $adv_sound_ratio $delta_t_between_analyses $div_damp_4th_order_switch $dissipative_heating $diff_h_smag_fac $shear_bg $damping_start_height_over_toa $damping_coeff_max
fi
cd - > /dev/null
