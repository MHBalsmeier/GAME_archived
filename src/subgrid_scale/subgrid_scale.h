/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

int calc_temp_diffusion_coeffs(State *, Config *, Irreversible_quantities *, Diagnostics *, double, Grid *);
int calc_mass_diffusion_coeffs(State *, Config *, Irreversible_quantities *, Diagnostics *, double, Grid *);
int hori_div_viscosity_eff(State *, Irreversible_quantities *, Grid *, Diagnostics *, Config *);
int hori_curl_viscosity_eff_rhombi(State *, Irreversible_quantities *, Grid *, Diagnostics *, Config *);
int hori_curl_viscosity_eff_triangles(State *, Irreversible_quantities *, Grid *, Dualgrid *dualgrid, Diagnostics *, Config *);
int vert_w_viscosity_eff(State *, Grid *, Diagnostics *, Irreversible_quantities *, double);
int vert_hor_mom_viscosity(State *, Irreversible_quantities *, Diagnostics *, Config *, Grid *, double);
int tke_update(Irreversible_quantities *, double);
