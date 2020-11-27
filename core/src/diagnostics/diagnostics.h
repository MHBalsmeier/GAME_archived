/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int calc_mass_diffusion_coeffs(State *, Config_info *, Scalar_field, Scalar_field);
int calc_temp_diffusion_coeffs(State *, Config_info *, Scalar_field, Scalar_field);
int calc_divv_term_viscosity_eff(State *, Config_info *, Scalar_field);
int calc_curl_term_viscosity_eff(State *, Config_info *, Scalar_field);
int remap_horpri2hordual_vector(Vector_field, int, int, double *, Grid *);
int vorticity_flux_horizontal(Vector_field, Curl_field, int, int, double *, Grid *);
int tangential_wind(Vector_field, int, int, double *, Grid *);
int remap_verpri2horpri_vector(Vector_field, int, int, double *, Grid *);
int vorticity_flux_vertical(Vector_field, Curl_field, int, int, double *, Grid *);
double spec_heat_cap_diagnostics_p(State *, int);
double spec_heat_cap_diagnostics_v(State *, int);
int pot_temp_diagnostics_dry(State *, Scalar_field);
double gas_constant_diagnostics(State *, int);
double calc_micro_density(double, double);
double calc_condensates_density_sum(int, Mass_densities);
int vertical_contravariant(Vector_field, int, int, Grid *, double *);
int horizontal_covariant(Vector_field, int, int, Grid *, double *);
int temperature_diagnostics_explicit(State *, State *, Diagnostics *, double);
double density_total(State *, int);
double density_gas(State *, int);
