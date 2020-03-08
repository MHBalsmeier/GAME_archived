imit_time=0
run_cfg_filename=test_0
run_cfg=run_configs/$run_cfg_filename
output_folder=output/$run_cfg_filename
if [ -d $output_folder ]
then
rm -r $output_folder
fi
mkdir $output_folder
./game $run_cfg $init_time
