/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

int forward_tendencies(State *, State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Interpolate_info *, Diffusion_info *, Config_info *, int);
int backward_tendencies(State *, State *, Interpolate_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);
int integrate_continuity_dry(State *, State *, Interpolate_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);
int integrate_tracers(State *, State *, Interpolate_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);
int integrate_temp_gas(State *, State *, Interpolate_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);
int integrate_entropy_density_gas(State *, State *, Interpolate_info *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Diffusion_info *, Config_info *, int);
