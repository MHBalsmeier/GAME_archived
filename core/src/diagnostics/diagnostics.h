/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

int calc_temp_diffusion_coeffs(State *, Config_info *, Irreversible_quantities *, Diagnostics *, double, Grid *);
int calc_mass_diffusion_coeffs(State *, Config_info *, Scalar_field, Scalar_field);
int hori_viscosity_eff(State *, Vector_field, Grid *, Diagnostics *, Config_info *, double);
int ver_viscosity_eff(State *, Vector_field, Grid *);
int remap_horpri2hordual_vector(Vector_field, int, int, double *, Grid *);
int vorticity_flux_horizontal_traditional(Vector_field, Curl_field, int, int, double *, Grid *);
int vorticity_flux_vertical(Vector_field, Curl_field, int, int, double *, Grid *, Dualgrid *);
int tangential_wind(Vector_field, int, int, double *, Grid *);
int remap_verpri2horpri_vector(Vector_field, int, int, double *, Grid *);
int pot_temp_diagnostics_dry(State *, Scalar_field);
double spec_heat_cap_diagnostics_p(State *, int, Config_info *);
double spec_heat_cap_diagnostics_v(State *, int, Config_info *);
double gas_constant_diagnostics(State *, int, Config_info *);
double calc_micro_density(double, double);
double calc_condensates_density_sum(int, Mass_densities);
int vertical_contravariant_corr(Vector_field, int, int, Grid *, double *);
int horizontal_covariant(Vector_field, int, int, Grid *, double *);
int temperature_step(State *, State *, State *, Diagnostics *, Config_info *, double, int);
double density_total(State *, int);
double density_gas(State *, int);
int entropy_density_step(State *, State *, State *, Diagnostics *, Config_info *, double);
int calc_uv_at_edge(Vector_field, Vector_field, Vector_field, Grid *);
