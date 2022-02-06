/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

#define DENSITY_WATER 1024.0
#define GRAVITY 9.8
#define KARMAN 0.4

double sink_velocity(int, double, double);
int calc_h2otracers_source_rates(double [], double [], double [], double [], double [], int, double, int, int [], int,
double [], double [], double [], int, double [], double [], double [], double [], int);
double water_vapour_density_from_rel_humidity(double, double, double);
double saturation_pressure_over_water(double);
double saturation_pressure_over_ice(double);
double c_p_cond(int, int, double);
double c_v_cond(int, int, double);
double phase_trans_heat(int, double);
double rel_humidity(double, double);
double entropy_constants_gas_lookup(int);
double mean_particle_masses_gas_lookup(int);
double spec_heat_capacities_v_gas_lookup(int);
double spec_heat_capacities_p_gas_lookup(int);
double specific_gas_constants_lookup(int);
double solve_specific_entropy_for_density(double, double);
double spec_entropy_from_temp(double, double);
double calc_o3_vmr(double);
