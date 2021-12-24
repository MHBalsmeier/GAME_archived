/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

int manage_rkhevi(State *, State *, Soil *, Grid *, Dualgrid *, State *, Diagnostics *, Forcings *, Irreversible_quantities *, Config *, double, double, int);
int moisturizer(State *, double, Diagnostics *, Irreversible_quantities *, Config *, Grid *, Soil *);
int vector_tendencies_expl(State *, State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Irreversible_quantities *, Config *, int, int, double);
int scalar_tendencies_expl(State *, State *, State *, Soil *, Grid *, double, Diagnostics *, Forcings *, Irreversible_quantities *, Config *, int, int);
int three_band_solver_gen_densitites(State *, State *, State *, Diagnostics *, Config *, double, Grid *);
int three_band_solver_ver_waves(State *, State *, State *, Diagnostics *, Config *, double, Grid *, int);
int manage_pressure_gradient(State *, Grid *, Dualgrid *, Diagnostics *, Forcings *,  Irreversible_quantities *, Config *);
