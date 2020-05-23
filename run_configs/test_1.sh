#!/bin/bash
operator=boss
run_span=400
write_out_interval=200
grid_props_file=/home/max/compiled/game/grids/B4L6T30000_M2_O1_OL4.nc
init_state_filename=test_1_B4L6T30000_M2_O1_OL4.grb2
init_state_file=/home/max/compiled/game/input/$init_state_filename
output_dir_base=/home/max/compiled/game/output
cfl_margin=0.2
dissipation=1
rad_on=1
add_comps_on=1
write_out_dry_mass_integral=0
write_out_entropy_integral=0
write_out_kinetic_energy_integral=0
write_out_potential_energy_integral=0
write_out_internal_energy_integral=0
source run.sh
