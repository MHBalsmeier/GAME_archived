/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

// The number of pressure_leveltic pressure levels for the pressure_leveltic output.
#define NO_OF_PRESSURE_LEVELS 6

// The number of flight levels for the flight_level output.
#define NO_OF_FLIGHT_LEVELS 7

int get_damping_layer_properties(double *, double *);
double get_impl_vert_compression_weight();
double get_expl_pgrad_weight();
double get_impl_w_vadv_weight();
double get_impl_u_vadv_weight();
double get_impl_t_adv_weight();
int ask_for_divergence_output(int *);
double entropy_constants_gas(int);
double mean_particle_masses_gas(int);
double spec_heat_capacities_v_gas(int);
double spec_heat_capacities_p_gas(int);
double specific_gas_constants(int);
int get_flight_levels(double []);
int get_pressure_levels(double []);
