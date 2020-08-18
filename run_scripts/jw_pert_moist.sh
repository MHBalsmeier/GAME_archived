#!/bin/bash
operator=MHB
overwrite_run_id=1
run_id=jw_pert_moist
run_span=43200
write_out_interval=900
grid_props_file=grids/B6L26T30000_O2_OL17_SCVT.nc
init_state_filename=test_5_B6L26T30000_O2_OL17_SCVT.nc
init_state_file=input/$init_state_filename
output_dir_base=output
cfl_margin=0.0
temperature_diff_h=1
temperature_diff_v=1
momentum_diff_h=1
momentum_diff_v=1
tracers_on=1
rad_on=0
radiation_delta_t=3600
write_out_mass_dry_integral=1
write_out_entropy_gas_integral=1
write_out_linearized_entropy_gas_integral=1
write_out_energy_integral=1
# relevant only for OMP
export OMP_NUM_THREADS=5
# relevant only for MPI
number_of_cpus=1
year=2000
month=1
day=1
hour=0
nwp_mode=0
# necessary only for data assimilation
ndvar_directory=/home/max/compiled/ndvar
source core/run.sh
