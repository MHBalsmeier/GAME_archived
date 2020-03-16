#!/bin/bash
if [ -d $output_folder ]
then
rm -r $output_folder
fi
mkdir $output_folder
./game $run_span $write_out_interval $geo_prop_file $init_state_file $output_folder
