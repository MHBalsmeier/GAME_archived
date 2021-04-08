/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

int manage_rkhevi(State *, State *, Extrapolation_info *, Grid *, Dualgrid *, Radiation *, State *, Diagnostics *, Forcings *, Irreversible_quantities *, Config_info *, double, double, int);
int vector_tendencies_expl(State *, State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Irreversible_quantities *, Config_info *, int, int, double);
int scalar_tendencies_expl(State *, State *, Grid *, Dualgrid *, double, Scalar_field, Diagnostics *, Forcings *, Irreversible_quantities *, Config_info *, int, Vector_field);
int three_band_solver_gen_densitites(State *, State *, State *, Diagnostics *, Config_info *, double, Grid *);
int three_band_solver_ver_waves(State *, State *, State *, Diagnostics *, Config_info *, double, Grid *, int);
int manage_pressure_gradient(State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Extrapolation_info *, Irreversible_quantities *, Config_info *);
