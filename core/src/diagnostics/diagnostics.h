/*
This source file is part of the Global Geophysical Modeling Frame (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int calc_temp_diffusion_coeffs(Scalar_field, Scalar_field, Scalar_field, Scalar_field);
int calc_mass_diffusion_coeffs(Scalar_field, Scalar_field, Scalar_field, Scalar_field);
int global_scalar_integrator(Scalar_field, Grid *, double *);
int recov_hor_par_curl(Curl_field, int, int, double *, Grid *);
int trsk_modified(Vector_field, Curl_field, int, int, double *, Grid *);
int recov_hor_par_pri(Vector_field, int, int, double *, Grid *);
int recov_hor_ver_pri(Vector_field, int, int, double *, Grid *);
int recov_ver_0_curl(Dual_vector_field, int, int, double *, Grid *);
int recov_ver_0_pri(Vector_field, int, int, double *, Grid *);
int recov_ver_1_curl(Dual_vector_field, int, int, double *, Grid *);
int recov_ver_1_pri(Vector_field, int, int, double *, Grid *);
double spec_heat_cap_diagnostics_p(double, double);
double spec_heat_cap_diagnostics_v(double, double);
int pot_temp_diagnostics(Scalar_field, Scalar_field, Tracer_densities, Scalar_field);
double pot_temp_diagnostics_single_value(double, double, double, double);
double exner_pressure_diagnostics_single_value(double, double, double);
int exner_pressure_diagnostics(Scalar_field, Scalar_field, Tracer_densities, Scalar_field);
double temperature_diagnostics_single_value(double, double);
int temperature_diagnostics(Scalar_field, Scalar_field, Tracer_densities, Scalar_field);
double gas_constant_diagnostics(double, double);
double calc_micro_density(double, double);
double calc_condensates_density_sum(int, int, Tracer_densities);
int vertical_contravariant_normalized(Vector_field, int, int, Grid *, double *);
int horizontal_covariant_normalized(Vector_field, int, int, Grid *, double *);
