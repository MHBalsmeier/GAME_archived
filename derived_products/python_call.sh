disp_level=0
disp_shortname=pres
fig_save_path=~/figs/game_output/$run_id
plot_integrals=1
if [ -d $fig_save_path ]
then
rm -r $fig_save_path
fi
mkdir $fig_save_path
echo plotting maps ...
python derived_products/create_movie_elements.py $run_span $write_out_interval $disp_level $disp_shortname $grid_props_file $fig_save_path $output_dir
echo done
if [ $plot_integrals -eq 1 ]
then
echo plotting integrals ...
python derived_products/create_integral_plots.py $fig_save_path $output_dir $write_out_dry_mass_integral $write_out_entropy_integral $write_out_energy_integral
echo done
fi
