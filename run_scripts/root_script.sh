#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

run_dir=../output/$run_id
if [ -d $run_dir ]
then
rm -r $run_dir
fi
mkdir $run_dir
cd $run_dir

if [ ! -f ../../build/game ]
then
echo "Executable game missing. Compile first. Aborting run."
cd - > /dev/null
exit 1
fi

if [ -f game ]
then
rm game
fi

cp ../../build/game .

if [ $valgrind_check -eq 0 ]
then
mpirun -np $number_of_cpus ./game $run_span $write_out_interval $cfl_margin $momentum_diff_h $momentum_diff_v $rad_on $operator $write_out_mass_integrals $write_out_rhotheta_integral $write_out_energy_integral $temperature_diff_h $radiation_delta_t $start_year $start_month $start_day $start_hour $temperature_diff_v $run_id $toa $orography_id $ideal_input_id $grib_output_switch $netcdf_output_switch $pressure_level_output_switch $model_level_output_switch $surface_output_switch $orography_layers $type_of_vertical_grid $assume_lte $adv_sound_ratio $delta_t_between_analyses $dissipative_heating $diff_h_smag_fac $shear_bg $damping_start_height_over_toa $damping_coeff_max $explicit_boundary_layer
fi
if [ $valgrind_check -eq 1 ]
then
valgrind ./game $run_span $write_out_interval $cfl_margin $momentum_diff_h $momentum_diff_v $rad_on $operator $write_out_mass_integrals $write_out_rhotheta_integral $write_out_energy_integral $temperature_diff_h $radiation_delta_t $start_year $start_month $start_day $start_hour $temperature_diff_v $run_id $toa $orography_id $ideal_input_id $grib_output_switch $netcdf_output_switch $pressure_level_output_switch $model_level_output_switch $surface_output_switch $orography_layers $type_of_vertical_grid $assume_lte $adv_sound_ratio $delta_t_between_analyses $dissipative_heating $diff_h_smag_fac $shear_bg $damping_start_height_over_toa $damping_coeff_max $explicit_boundary_layer
fi
cd - > /dev/null
