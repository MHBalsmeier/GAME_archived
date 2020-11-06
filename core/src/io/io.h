/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int set_grid_properties(Grid *, Dualgrid *, char[]);
int calc_delta_t(double, double *, Grid *);
int set_init_data(char[], State *);
int write_out(State *, double [], int, double, double, Diagnostics *, Forcings *, Grid *, Dualgrid *, char [], Io_config *);
int write_out_integral(State *, int, char [], Grid *, Dualgrid *, Diagnostics *, int);
int interpolation_t(State *, State *, State *, double, double, double);
