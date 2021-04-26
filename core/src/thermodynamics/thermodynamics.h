/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

int pot_temp_diagnostics_dry(State *, Scalar_field);
double spec_heat_cap_diagnostics_p(State *, int, Config_info *);
double spec_heat_cap_diagnostics_v(State *, int, Config_info *);
double gas_constant_diagnostics(State *, int, Config_info *);
double calc_micro_density(double, double);
double calc_condensates_density_sum(int, Mass_densities);
double density_total(State *, int);
double density_gas(State *, int);
int calc_diffusion_coeff(double, double, double, double, double *);
