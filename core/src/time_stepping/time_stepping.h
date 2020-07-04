/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int manage_time_stepping(State *, State *, double, Grid *, Dualgrid *, int, int, int, int, Scalar_field, double [], double [], State *, Vector_field, Scalar_field, Scalar_field, Vector_field, Scalar_field, Scalar_field, Vector_field, Vector_field, Scalar_field, Scalar_field, Curl_field, Curl_field, Vector_field, Vector_field, Vector_field, Vector_field, Scalar_field, Scalar_field, Scalar_field, Vector_field, Scalar_field, Scalar_field, Vector_field, Scalar_field, Vector_field, Scalar_field, Vector_field, Vector_field, Scalar_field, Scalar_field, Vector_field, Scalar_field, Vector_field);
int solve_lower_boundary(State *, Grid *);
int three_band_solver_hor(State *, State *, State *, double, Grid *);
int three_band_solver_ver_vel(State *, State *, State *, double, Grid *);
int three_band_solver_ver_den_dry(State *, State *, State *, double, Grid *);
int three_band_solver_ver_entropy_gas(State *, State *, State *, double, Grid *);
int three_band_solver_ver_tracers(State *, State *, State *, double, Grid *);
int vertical_sound_wave_solver(State *, State *, State *, double, Grid *, Vector_field, int);
