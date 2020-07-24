#!/bin/bash
operator=MHB
overwrite_run_id=1
source run_configs/.sh/determine_latest_synop_time.sh
run_id=$year$month$day$hour"_nwp"
echo $run_id
run_span=259200
write_out_interval=3600
grid_props_file=grids/B5L26T30000_O3_OL17_SCVT.nc
init_state_filename=$year$month$day$hour"_nwp_B5L26T30000_O3_OL17_SCVT.nc"
init_state_file=input/$init_state_filename
output_dir_base=output
cfl_margin=0.0
scalar_diffusion_on=1
momentum_diffusion_on=1
tracers_on=1
rad_on=1
radiation_delta_t=10800
write_out_mass_dry_integral=0
write_out_entropy_gas_integral=0
write_out_energy_integral=0
export OMP_NUM_THREADS=4
number_of_cpus=1
source core/run.sh
