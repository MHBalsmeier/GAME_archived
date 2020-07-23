/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int manage_time_stepping(State *, State *, Interpolate_info *, Grid *, Dualgrid *, Scalar_field, State *, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, double);
int solve_lower_boundary(State *, Grid *);
int three_band_solver_hor(State *, State *, State *, double, Grid *);
int three_band_solver_ver_sound_waves(State *, State *, State *, Diagnostics *, double, Grid *);
int three_band_solver_ver_vel_adv(State *, State *, State *, double, Grid *);
int three_band_solver_ver_den_dry(State *, State *, State *, double, Grid *);
int three_band_solver_ver_entropy_density_gas(State *, State *, State *, double, Grid *);
int three_band_solver_ver_tracers(State *, State *, State *, double, Grid *);
int setup_interpolation(State *, Interpolate_info *);
int update_interpolation(State *, State *, Interpolate_info *);
int set_interpolated_temperature(State *, State *, Interpolate_info *, int);
