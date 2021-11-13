/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

double spec_heat_cap_diagnostics_p(State *, int, Config *);
double spec_heat_cap_diagnostics_v(State *, int, Config *);
double gas_constant_diagnostics(State *, int, Config *);
double mean_particle_masses_gas(int);
double spec_heat_capacities_v_gas(int);
double spec_heat_capacities_p_gas(int);
double specific_gas_constants(int);
double calc_micro_density(double, double);
double calc_condensates_density_sum(int, Mass_densities);
double density_total(State *, int);
double density_gas(State *, int);
int calc_diffusion_coeff(double, double, double, double, double *);
int temperature_diagnostics(State *, Grid *, Diagnostics *);
