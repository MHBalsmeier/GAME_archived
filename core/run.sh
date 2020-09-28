#!/bin/bash
time_string=$(date --utc +%Y%m%d%H%M%S)
if [ $overwrite_run_id -eq 0 ]
then
run_id="$operator"_"$init_state_filename"_"$time_string"
fi
output_dir=$output_dir_base/$run_id
if [ -d $output_dir ]
then
rm -r $output_dir
fi
if [ ! -f $grid_props_file ]
then
echo "The file $grid_props_file does not exist. Please create and install the grid file before running the model."
echo "Aborting."
return
fi
if [ ! -f $init_state_file ]
then
echo "The file $init_state_file does not exist. Please create and install the initialization state before running the model."
echo "Aborting."
return
fi
mkdir $output_dir
if [ $nwp_mode -eq 1 ]
then
cd $ndvar_directory
source run.sh
cd - > /dev/null
fi
mpirun -np $number_of_cpus ./core/game $run_span $write_out_interval $grid_props_file $init_state_file $output_dir $cfl_margin $momentum_diff $rad_on $tracers_on $operator $write_out_mass_dry_integral $write_out_entropy_gas_integral $write_out_energy_integral $temperature_diff_h $radiation_delta_t $year $month $day $hour $temperature_diff_v $run_id $write_out_linearized_entropy_gas_integral
# valgrind ./core/game $run_span $write_out_interval $grid_props_file $init_state_file $output_dir $cfl_margin $momentum_diff $rad_on $tracers_on $operator $write_out_mass_dry_integral $write_out_entropy_gas_integral $write_out_energy_integral $temperature_diff_h $radiation_delta_t $year $month $day $hour $temperature_diff_v $run_id
