/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int forward_tendencies(State *, State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Interpolation_info *, Irreversible_quantities *, Config_info *, int);
int integrate_momentum(State *, State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Irreversible_quantities *, Config_info *, int);
int backward_tendencies(State *, Interpolation_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Irreversible_quantities *, Config_info *, int);
int integrate_generalized_densities(State *, Interpolation_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Irreversible_quantities *, Config_info *, int);
int three_band_solver_ver_hor_vel_adv(State *, State *, State *, double, Grid *);
int three_band_solver_gen_densitites(State *, State *, State *, Diagnostics *, double, Grid *);
int three_band_solver_ver_sound_waves(State *, Diagnostics *, Interpolation_info *, double, Grid *);
int manage_rkhevi(State *, State *, Interpolation_info *, Grid *, Dualgrid *, Scalar_field, State *, Diagnostics *, Forcings *, Irreversible_quantities *, Config_info *, double);
int setup_interpolation(State *, Interpolation_info *);
int update_interpolation(State *, State *, Interpolation_info *);
int manage_pressure_gradient(State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Interpolation_info *, Irreversible_quantities *, Config_info *);
