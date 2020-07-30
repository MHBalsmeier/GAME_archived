/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int calc_temp_diffusion_coeffs(Scalar_field, Scalar_field, Scalar_field, Scalar_field);
int calc_mass_diffusion_coeffs(Scalar_field, Scalar_field, Scalar_field, Scalar_field);
int global_scalar_integrator(Scalar_field, Grid *, double *);
int recov_hor_par_curl(Curl_field, int, int, double *, Grid *);
int recov_primal2dual(Vector_field, int, int, double *, Grid *);
int trsk_modified(Vector_field, Curl_field, int, int, double *, Grid *);
int recov_hor_par_pri(Vector_field, int, int, double *, Grid *);
int recov_hor_ver_pri(Vector_field, int, int, double *, Grid *);
int vertical_coriolis_gen(Vector_field, Curl_field, int, int, double *, Grid *);
double spec_heat_cap_diagnostics_p(double, double);
double spec_heat_cap_diagnostics_v(double, double);
double entropy_constant_diagnostics(double, double);
int pot_temp_diagnostics(State *, Scalar_field);
double gas_constant_diagnostics(double, double);
double calc_micro_density(double, double);
double calc_condensates_density_sum(int, Tracer_densities);
int vertical_contravariant_normalized(Vector_field, int, int, Grid *, double *);
int vertical_contravariant_normalized_h(Vector_field, int, int, Grid *, double *);
int horizontal_covariant_normalized(Vector_field, int, int, Grid *, double *);
int temperature_diagnostics(State *, State *);
