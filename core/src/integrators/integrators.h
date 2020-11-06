/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int forward_tendencies(State *, State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Interpolate_info *, Diffusion_info *, Config_info *, int);
int integrate_momentum(State *, State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);
int backward_tendencies(State *, State *, Interpolate_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);
int integrate_generalized_densities(State *, State *, Interpolate_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);
int three_band_solver_ver_vel_adv(State *, State *, State *, double, Grid *);
int three_band_solver_gen_densitites(State *, State *, State *, Diagnostics *, Config_info *, double, Grid *);
int three_band_solver_ver_sound_waves(State *, State *, State *, Diagnostics *, double, Grid *);
int manage_rkhevi(State *, State *, Interpolate_info *, Grid *, Dualgrid *, Scalar_field, State *, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, double);
int setup_interpolation(State *, Interpolate_info *);
int update_interpolation(State *, State *, Interpolate_info *);
int manage_pressure_gradient(State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Interpolate_info *, Diffusion_info *, Config_info *, int);
