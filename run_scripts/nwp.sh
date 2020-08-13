#!/bin/bash
operator=MHB
overwrite_run_id=1
source run_configs/.sh/determine_latest_synop_time.sh
run_id=$year$month$day$hour"_nwp"
run_span=0
write_out_interval=900
grid_props_file=grids/B6L26T30000_O3_OL17_SCVT.nc
init_state_filename=$year$month$day$hour"_nwp_B6L26T30000_O3_OL17_SCVT.nc"
init_state_file=input/$init_state_filename
output_dir_base=output
cfl_margin=0.0
temperature_diff_h=1
temperature_diff_v=1
momentum_diff_h=1
momentum_diff_v=1
tracers_on=1
rad_on=1
radiation_delta_t=10800
write_out_mass_dry_integral=0
write_out_entropy_gas_integral=0
write_out_energy_integral=0
# relevant only for OMP
export OMP_NUM_THREADS=5
# relevant only for MPI
number_of_cpus=1
nwp_mode=1
# necessary only for data assimilation
ndvar_directory=/home/max/compiled/ndvar
source core/run.sh
