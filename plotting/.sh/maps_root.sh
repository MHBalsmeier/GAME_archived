#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

if [ -d $fig_save_path ]
then
rm -r $fig_save_path
fi
mkdir $fig_save_path
echo plotting maps ...
number_of_variables=${#disp_shortname_list[@]}
for (( i=0; $i<$number_of_variables ; i=$(($i+1)) ))
do
python3 $game_home_dir/plotting/.py/plt_maps.py $run_span ${plot_intervals_list[$i]} ${disp_level_list[$i]} ${disp_shortname_list[$i]} $grid_props_file $fig_save_path $output_dir ${projections_list[$i]} $run_id ${uniform_colormap_list[$i]} ${scope_list[$i]} ${on_pressure_level_list[$i]} ${synoptical_time_mode[$i]} $init_year $init_month $init_day $init_hr
done
echo Finished plotting maps.
