#!/bin/bash
if [ -d $output_dir ]
then
rm -r $output_dir
fi
mkdir $output_dir
./core/game $run_span $write_out_interval $grid_props_file $init_state_file $output_dir $cfl_margin $dissipation $rad_on $add_comps_on
