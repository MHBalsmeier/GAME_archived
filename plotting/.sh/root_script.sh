if [ -d $fig_save_path ]
then
rm -r $fig_save_path
fi
mkdir $fig_save_path
echo plotting maps ...
python plotting/maps.py $run_span $write_out_interval $disp_level $disp_shortname $grid_props_file $fig_save_path $output_dir $projection
echo Finished plotting maps.
if [ $plot_integrals -eq 1 ]
then
echo plotting integrals ...
python plotting/integrals.py $fig_save_path $output_dir $write_out_mass_dry_integral $write_out_entropy_gas_integral $write_out_energy_integral
echo Finished plotting integrals.
fi
