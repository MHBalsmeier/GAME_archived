#!/bin/bash

run_id=ullrich_dry_rev
output_dir=/home/max/compiled/game_dev/output/$run_id
fig_save_path=/home/max/figs/game_output/$run_id
plot_mass_dry_integral=1
plot_entropy_gas_integral=1
plot_energy_integral=1
plot_linearized_entropy_gas_integral=1
source plotting/.sh/integrals_root.sh
