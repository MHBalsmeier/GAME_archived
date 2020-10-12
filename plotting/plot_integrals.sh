#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

game_home_dir=/home/max/compiled/game_dev
run_id=ullrich_dry_rev
output_dir=$game_home_dir/output/$run_id
fig_save_path=/home/max/figs/game_output/$run_id
plot_mass_dry_integral=1
plot_entropy_gas_integral=1
plot_energy_integral=1
plot_linearized_entropy_gas_integral=1
source $game_home_dir/plotting/.sh/integrals_root.sh
