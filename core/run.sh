#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

time_string=$(date --utc +%Y%m%d%H%M%S)

# this is the NWP case
if [ $ideal_input_id -eq "-1" ]
then
run_id=$start_year$start_month$start_day$start_hour"_nwp_B"$res_id"L"$number_of_layers"T"$toa"_O"$orography_id"_OL"$orography_layers"_SCVT"
else
orography_id=-1
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
mpirun -np $number_of_cpus $game_home_dir/core/game $run_span $write_out_interval $res_id $number_of_layers $cfl_margin $momentum_diff $rad_on $tracers_on $operator $write_out_mass_dry_integral $write_out_entropy_gas_integral $write_out_energy_integral $temperature_diff_h $radiation_delta_t $start_year $start_month $start_day $start_hour $temperature_diff_v $run_id $write_out_linearized_entropy_gas_integral $toa $orography_layers $orography_id $ideal_input_id $mass_dry_diff_h $mass_dry_diff_v $grib_output_switch $netcdf_output_switch $synop_output_mode $aviation_output_mode
fi
if [ $valgrind_check -eq 1 ]
then
valgrind $game_home_dir/core/game $run_span $write_out_interval $res_id $number_of_layers $cfl_margin $momentum_diff $rad_on $tracers_on $operator $write_out_mass_dry_integral $write_out_entropy_gas_integral $write_out_energy_integral $temperature_diff_h $radiation_delta_t $start_year $start_month $start_day $start_hour $temperature_diff_v $run_id $write_out_linearized_entropy_gas_integral $toa $orography_layers $orography_id $ideal_input_id $mass_dry_diff_h $mass_dry_diff_v $grib_output_switch $netcdf_output_switch $synop_output_mode $aviation_output_mode
fi
cd - > /dev/null
