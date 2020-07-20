/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int grad(Scalar_field, Vector_field, Grid *);
int calc_pot_vort(Vector_field, Scalar_field, Curl_field, Grid *, Dualgrid *);
int calc_rel_vort(Vector_field, Curl_field, Grid *, Dualgrid *);
int coriolis_gen(Vector_field, Dual_vector_field, Vector_field, Grid *);
int kinetic_energy(Vector_field, Scalar_field, Grid *);
int divv(Vector_field, Scalar_field, Grid *, int);
int divv_h(Vector_field, Scalar_field, Grid *);
int explicit_momentum_tendencies(State *, State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Interpolate_info *, Diffusion_info *, Config_info *, int);
int calc_partially_implicit_divvs(State *, State *, Interpolate_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);
int scalar_times_scalar(Scalar_field, Scalar_field, Scalar_field);
int scalar_times_vector(Scalar_field, Vector_field, Vector_field, Grid *);
int scalar_times_vector_scalar_h_v(Scalar_field, Scalar_field, Vector_field, Vector_field, Grid *);
int scalar_times_vector_vector_h_v(Scalar_field, Vector_field, Vector_field, Vector_field, Grid *);
int scalar_times_vector_v_column(double [], double [], double []);
int linear_combine_two_states(State *, State *, State *, double, double);
int dissipation(Vector_field, Scalar_field, Vector_field, Scalar_field, Grid *);
int divv_v_columns(double [], double [], int, Grid *);
int grad_v_vector_column_to_vector_points(double [], double [], int, Grid *);
int grad_v_scalar_column(double [], double [], int, Grid *);
int grad_v_vector_column_to_scalar_points(double [], double [], int, Grid *);
int grad_v_scalar_column_to_scalar_points(double [], double [], int, Grid *);
