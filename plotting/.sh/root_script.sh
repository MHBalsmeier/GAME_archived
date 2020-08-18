if [ -d $fig_save_path ]
then
rm -r $fig_save_path
fi
mkdir $fig_save_path
echo plotting maps ...
number_of_variables=${#disp_shortname_list[@]}
for (( i=0; $i<$number_of_variables ; i=$(($i+1)) ))
do
python plotting/maps.py $run_span ${plot_intervals_list[$i]} ${disp_level_list[$i]} ${disp_shortname_list[$i]} $grid_props_file $fig_save_path $output_dir ${projections_list[$i]} $run_id ${uniform_colormap_list[$i]} ${scope_list[$i]}
done
echo Finished plotting maps.
if [ $plot_integrals -eq 1 ]
then
echo plotting integrals ...
python plotting/integrals.py $fig_save_path $output_dir $write_out_mass_dry_integral $write_out_entropy_gas_integral $write_out_energy_integral $run_id $write_out_linearized_entropy_gas_integral
echo Finished plotting integrals.
fi
