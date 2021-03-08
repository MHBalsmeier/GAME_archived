/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

// The number of pressure levels for the pressure level output.
#define NO_OF_PRESSURE_LEVELS 6

double get_impl_thermo_weight();
double get_t_vadv_parameter();
int ask_for_divergence_output(int *);
double entropy_constants_gas(int);
double mean_particle_masses_gas(int);
double spec_heat_capacities_v_gas(int);
double spec_heat_capacities_p_gas(int);
double specific_gas_constants(int);
int get_pressure_levels(double []);
