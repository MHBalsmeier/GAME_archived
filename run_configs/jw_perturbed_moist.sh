#!/bin/bash
operator=Boss
overwrite_run_id=1
run_id=jw_perturbed_moist
run_span=900
write_out_interval=300
grid_props_file=/home/max/compiled/game/grids/B4L6T30000_M2_O2_OL4.nc
init_state_filename=test_4_B4L6T30000_M2_O2_OL4.grb2
init_state_file=/home/max/compiled/game/input/$init_state_filename
output_dir_base=/home/max/compiled/game/output
cfl_margin=0.2
dissipation=1
rad_on=0
add_comps_on=1
write_out_dry_mass_integral=0
write_out_entropy_integral=0
write_out_energy_integral=0
source run.sh
