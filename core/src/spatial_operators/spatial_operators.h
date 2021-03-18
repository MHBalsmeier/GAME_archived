/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

int grad_hor_cov(Scalar_field, Vector_field, Grid *);
int grad_vert_cov(Scalar_field, Vector_field, Grid *);
int grad_oro_corr(Vector_field, Grid *);
int grad_cov(Scalar_field, Vector_field, Grid *);
int grad(Scalar_field, Vector_field, Grid *);
int grad_hor(Scalar_field, Vector_field, Grid *);
int calc_pot_vort(Vector_field, Scalar_field, Diagnostics *, Grid *, Dualgrid *);
int add_f_to_rel_vort(Curl_field, Curl_field, Dualgrid *);
int calc_rel_vort(Vector_field, Curl_field, Grid *, Dualgrid *);
int vorticity_flux(Vector_field, Dual_vector_field, Vector_field, Grid *, Dualgrid *);
int kinetic_energy(Vector_field, Scalar_field, Grid *);
int divv_h(Vector_field, Scalar_field, Grid *);
int add_vertical_divv(Vector_field, Scalar_field, Grid *);
int scalar_times_scalar(Scalar_field, Scalar_field, Scalar_field);
int scalar_times_vector(Scalar_field, Vector_field, Vector_field, Grid *);
int scalar_times_vector_scalar_h_v(Scalar_field, Scalar_field, Vector_field, Vector_field, Grid *);
int vector_times_vector(Vector_field, Vector_field, Vector_field);
int linear_combine_two_states(State *, State *, State *, double, double);
int momentum_diff_diss(State *, Diagnostics *, Irreversible_quantities*, Config_info *, Grid *, Dualgrid *, double);
int inner_product(Vector_field, Vector_field, Scalar_field, Grid *, int);
int curl_of_vorticity(Curl_field, Vector_field, Grid *, Dualgrid *, Config_info *);
int calc_horizontal_shear(State *, Diagnostics *, Grid *);
