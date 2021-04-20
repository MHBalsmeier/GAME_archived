/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

int calc_temp_diffusion_coeffs(State *, Config_info *, Irreversible_quantities *, Diagnostics *, double, Grid *);
int hori_div_viscosity_eff(State *, Irreversible_quantities *, Grid *, Diagnostics *, Config_info *, double);
int hori_curl_viscosity_eff(State *, Irreversible_quantities *, Grid *, Diagnostics *, Config_info *, double);
int vert_w_viscosity_eff(State *, Grid *, Diagnostics *, double);
int remap_verpri2horpri_vector(Vector_field, int, int, double *, Grid *);
int pot_temp_diagnostics_dry(State *, Scalar_field);
double spec_heat_cap_diagnostics_p(State *, int, Config_info *);
double spec_heat_cap_diagnostics_v(State *, int, Config_info *);
double gas_constant_diagnostics(State *, int, Config_info *);
double calc_micro_density(double, double);
double calc_condensates_density_sum(int, Mass_densities);
double density_total(State *, int);
double density_gas(State *, int);
int calc_diffusion_coeff(double, double, double, double, double *);
