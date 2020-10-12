/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int manage_time_stepping(State *, State *, Interpolate_info *, Grid *, Dualgrid *, Scalar_field, State *, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, double);
int setup_interpolation(State *, Interpolate_info *);
int update_interpolation(State *, State *, Interpolate_info *);
int manage_pressure_gradient(State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Interpolate_info *, Diffusion_info *, Config_info *, int);
