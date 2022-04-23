/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

int hori_div_viscosity(State *, Irreversible_quantities *, Grid *, Diagnostics *, Config *);
int hori_curl_viscosity_rhombi(State *, Irreversible_quantities *, Grid *, Diagnostics *, Config *);
int hori_curl_viscosity_triangles(State *, Irreversible_quantities *, Grid *, Dualgrid *dualgrid, Diagnostics *, Config *);
int vert_hori_mom_viscosity(State *, Irreversible_quantities *, Diagnostics *, Config *, Grid *, double);
int vert_vert_mom_viscosity(State *, Grid *, Diagnostics *, Irreversible_quantities *, double);
int temp_diffusion_coeffs(State *, Config *, Irreversible_quantities *, Diagnostics *, double, Grid *);
int mass_diffusion_coeffs(State *, Config *, Irreversible_quantities *, Diagnostics *, double, Grid *);
int tke_update(Irreversible_quantities *, double, State *, Diagnostics *, Grid *);
int update_sfc_turb_quantities(State *, Grid *, Diagnostics *, Config *, double);
int pbl_wind_tendency(State *, Diagnostics *, Irreversible_quantities *, Grid *, double);
double momentum_flux_resistance(double, double, double, double, double);
