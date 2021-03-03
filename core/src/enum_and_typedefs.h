/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include <math.h>
// some fundamental constants
#define K_B (1.380649e-23)
#define P_0 100000.0
#define SECONDS_PER_HOUR 3600
#define RHO_WATER 1024.0
#define OMEGA (7.292115e-5)
#define RADIUS 6371000.789927
#define FOOT 0.3048
#define EPSILON_SECURITY (1e-10)

enum grid_integers {
// This determines the horizontal resolution.
RES_ID = 5,
// This has to conform with the grid file and the initialization state file.
NO_OF_LAYERS = 26,
// The number of layers affected by orography. This also has to conform with the grid file and the initialization state file.
NO_OF_GASEOUS_CONSTITUENTS = 1,
NO_OF_CONDENSED_CONSTITUENTS = 0,
// Nothing should be changed by the user below this line.
NO_OF_CONSTITUENTS = (NO_OF_GASEOUS_CONSTITUENTS + NO_OF_CONDENSED_CONSTITUENTS),
NO_OF_BASIC_TRIANGLES = 20,
NO_OF_PENTAGONS = 12,
NO_OF_HEXAGONS = (int) (10*(pow(2, 2*RES_ID) - 1)),
NO_OF_EDGES = 3*NO_OF_BASIC_TRIANGLES/2,
NO_OF_LEVELS = NO_OF_LAYERS + 1,
NO_OF_SCALARS_H = NO_OF_PENTAGONS + NO_OF_HEXAGONS,
NO_OF_VECTORS_H = (5*NO_OF_PENTAGONS/2 + 6/2*NO_OF_HEXAGONS),
NO_OF_H_VECTORS = NO_OF_LAYERS*NO_OF_VECTORS_H,
NO_OF_V_VECTORS = NO_OF_LEVELS*NO_OF_SCALARS_H,
NO_OF_VECTORS_PER_LAYER = NO_OF_VECTORS_H + NO_OF_SCALARS_H,
NO_OF_TRIANGLES = (int) (NO_OF_BASIC_TRIANGLES*(pow(4, RES_ID))),
NO_OF_SCALARS = NO_OF_SCALARS_H*NO_OF_LAYERS,
NO_OF_VECTORS = NO_OF_H_VECTORS + NO_OF_V_VECTORS,
NO_OF_DUAL_SCALARS_H = NO_OF_TRIANGLES,
NO_OF_DUAL_H_VECTORS = NO_OF_LEVELS*NO_OF_VECTORS_H,
NO_OF_DUAL_V_VECTORS = NO_OF_LAYERS*NO_OF_DUAL_SCALARS_H,
NO_OF_DUAL_VECTORS_PER_LAYER = NO_OF_VECTORS_H + NO_OF_DUAL_SCALARS_H,
NO_OF_DUAL_SCALARS = NO_OF_LEVELS*NO_OF_DUAL_SCALARS_H,
NO_OF_DUAL_VECTORS = NO_OF_DUAL_H_VECTORS + NO_OF_DUAL_V_VECTORS};

typedef double Scalar_field[NO_OF_SCALARS];
typedef double Vector_field[NO_OF_VECTORS];
typedef double Dual_vector_field[NO_OF_DUAL_VECTORS];
typedef double Curl_field[NO_OF_LAYERS*2*NO_OF_VECTORS_H + NO_OF_VECTORS_H];
// all constituents have a mass density
typedef double Mass_densities[NO_OF_CONSTITUENTS*NO_OF_SCALARS];
// only the gaseous constituents have an entropy density
typedef double Entropy_densities[NO_OF_GASEOUS_CONSTITUENTS*NO_OF_SCALARS];
// only the condensed constituents have a density x temperature field
typedef double Condensed_density_temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS];

// Contains properties of the primal grid.
typedef struct grid {
int no_of_oro_layers;
Vector_field normal_distance;
Scalar_field volume;
Vector_field area;
Scalar_field z_scalar;
Vector_field z_vector;
Scalar_field gravity_potential;
Vector_field gravity_m;
Vector_field slope;
int trsk_indices[10*NO_OF_VECTORS_H];
int trsk_modified_curl_indices[10*NO_OF_VECTORS_H];
int from_index[NO_OF_VECTORS_H];
int to_index[NO_OF_VECTORS_H];
int adjacent_vector_indices_h[6*NO_OF_SCALARS_H];
int adjacent_signs_h[6*NO_OF_SCALARS_H];
int density_to_rhombus_indices[4*NO_OF_VECTORS_H];
int no_of_shaded_points_scalar[NO_OF_SCALARS_H];
int no_of_shaded_points_vector[NO_OF_VECTORS_H];
double latitude_scalar[NO_OF_SCALARS_H];
double longitude_scalar[NO_OF_SCALARS_H];
double volume_ratios[2*NO_OF_SCALARS];
double inner_product_weights[8*NO_OF_SCALARS];
double direction[NO_OF_VECTORS_H];
double density_to_rhombus_weights[4*NO_OF_VECTORS_H];
double trsk_weights[10*NO_OF_VECTORS_H];
double remap_horpri2hordual_vector_weights[2*NO_OF_DUAL_H_VECTORS];
double stretching_parameter;
double mean_area_edge;
} Grid;

// Contains properties of the dual grid.
typedef struct dualgrid {
Curl_field area;
Dual_vector_field normal_distance;
int from_index[NO_OF_VECTORS_H];
int to_index[NO_OF_VECTORS_H];
int adjacent_vector_indices_h[3*NO_OF_DUAL_SCALARS_H];
int vorticity_indices[4*NO_OF_VECTORS_H];
int vorticity_signs[4*NO_OF_VECTORS_H];
double f_vec[3*NO_OF_VECTORS_H];
} Dualgrid;

typedef struct state {
Scalar_field temperature_gas;
Vector_field velocity_gas;
// density order: solid, liquid, vapour
Mass_densities mass_densities;
Entropy_densities entropy_densities;
Condensed_density_temperatures condensed_density_temperatures;
} State;

// Collects diagnostic quantities. Note: in fact, forcings are also diagnostic quantities.
typedef struct diagnostics {
Vector_field flux_density;
Scalar_field flux_density_divv;
Vector_field temperature_gradient;
Scalar_field pressure_gradient_1_prefactor;
Scalar_field temperature_gas_explicit;
Curl_field rel_vort;
Curl_field pot_vort;
Scalar_field c_g_v_field;
Scalar_field c_g_p_field;
Scalar_field e_kin;
Vector_field cpgradt;
Vector_field tgrads_component;
Vector_field tgrads;
Scalar_field velocity_gas_divv;
Vector_field curl_of_vorticity_m;
Scalar_field scalar_field_placeholder;
Vector_field velocity_gen;
} Diagnostics;

// Collects forcings.
typedef struct forcings {
Vector_field pressure_gradient_acc_expl;
Vector_field e_kin_grad;
Vector_field pot_vort_tend;
} Forcings;

// Info on the run configuration is collected here.
typedef struct config_info {
int totally_first_step_bool;
int mass_diff_h;
int mass_diff_v;
int temperature_diff_h;
int temperature_diff_v;
int momentum_diff_h;
int momentum_diff_v;
int rad_on;
int rad_update;
int assume_lte;
int adv_sound_ratio;
int nwp_mode;
int delta_t_between_analyses;
double div_damp_coeff;
} Config_info;

// This is necessary for stability of horizontally propagating sound waves.
typedef struct extrapolation_info {
Vector_field cpgradt_m_old;
} Extrapolation_info;

// Contains everything on turbulence parametrizations as well as constituent-related quantities.
typedef struct irreversible_quantities {
Scalar_field temperature_diffusion_heating;
Vector_field friction_acc;
Scalar_field heating_diss;
Scalar_field scalar_diffusion_coeff_numerical_h;
Scalar_field scalar_diffusion_coeff_numerical_v;
Vector_field mass_diffusion_flux_density;
Scalar_field mass_diffusion_source_rate;
Scalar_field pressure_gradient_decel_factor;
double constituent_mass_source_rates[NO_OF_CONSTITUENTS*NO_OF_SCALARS];
double constituent_heat_source_rates[NO_OF_CONSTITUENTS*NO_OF_SCALARS];
Scalar_field viscosity_eff;
Vector_field velocity_grad_div;
} Irreversible_quantities;

// Info on input and output is collected here.
typedef struct io_config {
int grib_output_switch;
int netcdf_output_switch;
int pressure_level_output_switch;
int flight_level_output_switch;
int model_level_output_switch;
int surface_output_switch;
} Io_config;




