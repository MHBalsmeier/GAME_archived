/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int grad(Scalar_field, Vector_field, Grid *);
int curl(Vector_field, Curl_field, Grid *, Dualgrid *);
int coriolis_gen(Vector_field, Dual_vector_field, Vector_field, Grid *);
int kinetic_energy(Vector_field, Scalar_field, Grid *);
int divv(Vector_field, Scalar_field, Grid *, int);
int divv_h(Vector_field, Scalar_field, Grid *);
int explicit_momentum_tendencies(State *, State *, Grid *, Dualgrid *, int, int, Scalar_field, Scalar_field, Vector_field, Vector_field, Scalar_field, Scalar_field, Curl_field, Curl_field, Vector_field, Vector_field, Vector_field, Vector_field, Scalar_field, Scalar_field, Scalar_field, Vector_field);
int calc_partially_implicit_divvs(State *, State *, State *, Grid *, Dualgrid *, int, int, int, double, int, Scalar_field, int, double [], double [], Vector_field, Scalar_field, Scalar_field, Vector_field, Scalar_field, Scalar_field, Vector_field, Scalar_field, Scalar_field, Scalar_field, Vector_field, Scalar_field, Vector_field, Scalar_field, Vector_field, Vector_field, Scalar_field, Scalar_field, Vector_field, Scalar_field);
int scalar_times_scalar(Scalar_field, Scalar_field, Scalar_field);
int scalar_times_vector(Scalar_field, Vector_field, Vector_field, Grid *);
int scalar_times_vector_scalar_h_v(Scalar_field, Scalar_field, Vector_field, Vector_field, Grid *);
int scalar_times_vector_vector_h_v(Scalar_field, Scalar_field, Vector_field, Vector_field, Grid *);
int linear_combine_two_states(State *, State *, State *, double, double);
int linear_combine_two_states_scalars(State *, State *, State *, double, double);
int dissipation(Vector_field, Scalar_field, Vector_field, Scalar_field, Grid *);
int set_state_to_zero(State *);
