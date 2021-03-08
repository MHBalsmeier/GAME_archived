/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

int set_grid_properties(Grid *, Dualgrid *, char[]);
int calc_delta_t_and_related(double, double *, Grid *, Dualgrid *, State *, Config_info *);
int set_init_data(char[], State *, Grid *);
int write_out(State *, double [], int, double, double, Diagnostics *, Forcings *, Grid *, Dualgrid *, char [], Io_config *, Config_info *);
int write_out_integral(State *, int, char [], Grid *, Dualgrid *, Diagnostics *, int);
int interpolation_t(State *, State *, State *, double, double, double);
int epv_diagnostics(Curl_field, Scalar_field, Scalar_field, Grid *, Dualgrid *);
int interpolate_to_ll(double [], double [], Grid *);
