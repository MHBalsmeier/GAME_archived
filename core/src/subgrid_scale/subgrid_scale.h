/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

int calc_temp_diffusion_coeffs(State *, Config_info *, Irreversible_quantities *, Diagnostics *, double, Grid *);
int hori_div_viscosity_eff(State *, Irreversible_quantities *, Grid *, Diagnostics *, Config_info *, double);
int hori_curl_viscosity_eff(State *, Irreversible_quantities *, Grid *, Diagnostics *, Config_info *, double);
int vert_w_viscosity_eff(State *, Grid *, Diagnostics *, double);
int vert_hor_mom_viscosity(State *, Irreversible_quantities *, Diagnostics *, Config_info *, Grid *, double);
