#!/bin/bash

time_string=$(date --utc +%Y%m%d%H%M%S)
if [ $overwrite_run_id -eq 0 ]
then
run_id="$operator"_"$init_state_filename"_"$time_string"
fi
output_dir=$game_home_dir/output/$run_id
if [ -d $output_dir ]
then
rm -r $output_dir
fi
mkdir $output_dir
if [ $nwp_mode -eq 1 ]
then
cd $ndvar_directory
source run.sh
cd - > /dev/null
fi
cd $game_home_dir
mpirun -np $number_of_cpus $game_home_dir/core/game $run_span $write_out_interval $res_id $number_of_layers $cfl_margin $momentum_diff $rad_on $tracers_on $operator $write_out_mass_dry_integral $write_out_entropy_gas_integral $write_out_energy_integral $temperature_diff_h $radiation_delta_t $start_year $start_month $start_day $start_hour $temperature_diff_v $run_id $write_out_linearized_entropy_gas_integral $toa $orography_layers $orography_id $test_id $mass_dry_diff_h $mass_dry_diff_v
cd - > /dev/null
