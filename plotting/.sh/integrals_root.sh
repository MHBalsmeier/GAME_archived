#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

echo plotting integrals ...
python3 $game_home_dir/plotting/.py/plt_integrals.py $fig_save_path $output_dir $plot_mass_dry_integral $plot_entropy_gas_integral $plot_energy_integral $run_id $plot_linearized_entropy_gas_integral
echo Finished plotting integrals.
