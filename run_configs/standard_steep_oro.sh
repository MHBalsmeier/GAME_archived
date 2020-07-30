#!/bin/bash
operator=MHB
overwrite_run_id=1
run_id=standard_steep_oro
run_span=43200
write_out_interval=900
grid_props_file=grids/B5L26T30000_O1_OL17_SCVT.nc
init_state_filename=test_1_B5L26T30000_O1_OL17_SCVT.nc
init_state_file=input/$init_state_filename
output_dir_base=output
cfl_margin=0.0
scalar_diffusion_on=0
momentum_diffusion_on=0
tracers_on=0
rad_on=0
radiation_delta_t=3600
write_out_mass_dry_integral=1
write_out_entropy_gas_integral=1
write_out_energy_integral=1
export OMP_NUM_THREADS=4
number_of_cpus=1
year=2000
month=1
day=1
hour=0
nwp_mode=0
# necessary only for data assimilation
ndvar_directory=/home/max/compiled/ndvar
source core/run.sh
