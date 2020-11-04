/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include <math.h>
// some fundamental constants
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define M_V 0.0180152
#define EPSILON (M_V/M_D)
// #define R (N_A*K_B)
#define R 8.314463
// #define R_D (R/M_D)
#define R_D 287.057811
// #define R_V (R/M_V)
#define R_V 461.524879
#define C_D_P 1005.0
// #define C_D_V (C_D_P - R_D)
#define C_D_V 717.942189
#define C_V_P 1858.0
// #define C_V_V (C_V_P - R_V)
#define C_V_V 1396.475121
// #define DELTA_C_V_P (C_V_P - C_D_P)
#define DELTA_C_V_P 853.000000
#define P_0 100000.0
#define SECONDS_PER_HOUR 3600
#define RHO_WATER 1024.0
#define OMEGA (7.292115e-5)
#define SEMIMAJOR 6378137.0
#define SEMIMINOR 6356752.314
// #define RADIUS pow(SEMIMAJOR*SEMIMAJOR*SEMIMINOR, 1.0/3.0)
#define RADIUS 6371000.789927
#define H_BAR (1.054571817e-34)
// #define ENTROPY_CONSTANT_D (0.4*C_D_P*log(K_B/P_0*pow(M_D/N_A*K_B*exp(5.0/3)/(2*M_PI*H_BAR*H_BAR), 1.5)))
#define ENTROPY_CONSTANT_D 1566.752670
// #define ENTROPY_CONSTANT_V (0.4*C_V_P*log(K_B/P_0*pow(M_V/N_A*K_B*exp(5.0/3)/(2*M_PI*H_BAR*H_BAR), 1.5)))
#define ENTROPY_CONSTANT_V 2367.178359
#define EPSILON_SECURITY (1e-10)
#define FOOT 0.3048

enum grid_integers {
// This determines the horizontal resolution.
RES_ID = 5,
// The number of layers. TOA is determined by the grid file. This has to conform with the grid file and the initialization state file
NO_OF_LAYERS = 26,
// The number of layers affected by orography. This also has to conform with the grid file and the initialization state file
NO_OF_ORO_LAYERS = 17,
NO_OF_TRACERS = 3,
NO_OF_CONDENSED_TRACERS = 2,
NO_OF_SOLID_TRACERS = 1,
// Nothing may be changed below this line. These are fundamentals properties of the grid.
NO_OF_GASEOUS_TRACERS = NO_OF_TRACERS - NO_OF_CONDENSED_TRACERS,
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
NO_OF_DUAL_VECTORS = NO_OF_DUAL_H_VECTORS + NO_OF_DUAL_V_VECTORS,
TRIANGLES_PER_FACE = NO_OF_TRIANGLES/NO_OF_BASIC_TRIANGLES,
POINTS_PER_EDGE = (int) (pow(2, RES_ID) - 1),
SCALAR_POINTS_PER_INNER_FACE = (int) (0.5*(pow(2, RES_ID) - 2)*(pow(2, RES_ID) - 1)),
VECTOR_POINTS_PER_INNER_FACE = (int) (1.5*(pow(2, RES_ID) - 1)*pow(2, RES_ID))};

typedef double Scalar_field[NO_OF_SCALARS];
typedef double Vector_field[NO_OF_VECTORS];
typedef double Dual_vector_field[NO_OF_DUAL_VECTORS];
typedef double Curl_field[NO_OF_LAYERS*2*NO_OF_VECTORS_H + NO_OF_VECTORS_H];
typedef double Tracer_densities[NO_OF_TRACERS*NO_OF_SCALARS];
typedef double Tracer_entropy_densities[NO_OF_GASEOUS_TRACERS*NO_OF_SCALARS];
typedef double Tracer_density_temperatures[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS];

// Contains properties of the primal grid.
typedef struct grid {
Vector_field normal_distance;
Scalar_field volume;
Vector_field area;
Scalar_field z_scalar;
Vector_field z_vector;
Scalar_field gravity_potential;
Vector_field gravity_m;
Vector_field slope;
int trsk_modified_velocity_indices[10*NO_OF_VECTORS_H];
int trsk_modified_curl_indices[10*NO_OF_VECTORS_H];
int from_index[NO_OF_VECTORS_H];
int to_index[NO_OF_VECTORS_H];
int adjacent_vector_indices_h[6*NO_OF_SCALARS_H];
int adjacent_signs_h[6*NO_OF_SCALARS_H];
int density_to_rhombus_indices[4*NO_OF_VECTORS_H];
double volume_ratios[2*NO_OF_SCALARS];
double e_kin_weights[8*NO_OF_SCALARS];
double direction[NO_OF_VECTORS_H];
double density_to_rhombus_weights[4*NO_OF_VECTORS_H];
double trsk_modified_weights[10*NO_OF_VECTORS_H];
double recov_ver_weight[6*NO_OF_LEVELS*NO_OF_SCALARS_H];
double recov_primal2dual_weights[2*NO_OF_DUAL_H_VECTORS];
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
int h_curl_indices[4*NO_OF_VECTORS_H];
int h_curl_signs[4*NO_OF_VECTORS_H];
double f_vec[3*NO_OF_VECTORS_H];
} Dualgrid;

typedef struct state {
Scalar_field density_dry;
Scalar_field temperature_gas;
Scalar_field entropy_density_dry;
Vector_field velocity_gas;
// density order: solid, liquid, vapour
Tracer_densities tracer_densities;
Tracer_entropy_densities tracer_entropy_densities;
Tracer_density_temperatures tracer_density_temperatures;
} State;

// Collects diagnostic quantities. Note: in fact, forcings are also diagnostic quantities.
typedef struct diagnostics {
Vector_field mass_dry_flux_density;
Vector_field temperature_gradient;
Scalar_field specific_entropy_dry;
Scalar_field specific_entropy_vapour;
Scalar_field pressure_gradient_1_dry_prefactor;
Scalar_field pressure_gradient_1_vapour_prefactor;
Scalar_field temperature_gas_explicit;
Curl_field rel_vort;
Curl_field pot_vort;
Scalar_field c_h_v_field;
Scalar_field c_h_p_field;
Scalar_field e_kin_h;
// nabla h
Vector_field pressure_gradient_0_m;
// temperature times nabla s
Vector_field pressure_gradient_1_dry;
Vector_field pressure_gradient_1_vapour;
Vector_field entropy_dry_flux_density;
// This is for the momentum diffusion.
Scalar_field velocity_gas_divv;
Vector_field curl_of_vorticity_m;
} Diagnostics;

// Collects forcings.
typedef struct forcings {
Scalar_field mass_dry_flux_density_divv;
Scalar_field entropy_dry_flux_density_divv;
Scalar_field temperature_gas_flux_divv_h;
Vector_field pressure_gradient_acc;
Vector_field e_kin_h_grad;
Vector_field pot_vort_tend;
} Forcings;

// Info on the run configuration is collected here.
typedef struct config_info {
int totally_first_step_bool;
int mass_dry_diff_h;
int mass_dry_diff_v;
int temperature_diff_h;
int temperature_diff_v;
int momentum_diff;
int tracers_on;
int phase_transitions_on;
int rad_on;
int rad_update;
} Config_info;

// This is necessary for ensuring cancellation of energetically important terms, see Gassmann and Herzog.
typedef struct interpolate_info {
Vector_field pressure_gradient_0_old_m;
Vector_field pressure_gradient_1_old;
} Interpolate_info;

// Contains everything on turbulence parametrizations as well as tracer-related quantities.
typedef struct diffusion_info {
Vector_field temperature_diffusive_flux_density;
Scalar_field temperature_diffusion_heating;
Vector_field friction_acc;
Scalar_field heating_diss;
Scalar_field scalar_diffusion_coeff_numerical_h;
Scalar_field scalar_diffusion_coeff_numerical_v;
Vector_field mass_dry_diffusion_flux_density;
Scalar_field mass_dry_diffusion_source_rate;
Scalar_field tracer_density;
Scalar_field tracer_entropy_density;
Vector_field tracer_velocity;
Vector_field tracer_flux_density;
Scalar_field tracer_flux_density_divv;
Scalar_field tracer_density_temperature;
Vector_field tracer_temperature_flux_density;
Scalar_field tracer_temperature_flux_density_divv;
Scalar_field pressure_gradient_decel_factor;
double tracer_mass_source_rates[NO_OF_TRACERS*NO_OF_SCALARS];
double tracer_heat_source_rates[NO_OF_TRACERS*NO_OF_SCALARS];
Scalar_field divv_term_viscosity_eff;
Scalar_field curl_term_viscosity_eff;
} Diffusion_info;

// Info on input and output is collected here.
typedef struct io_config {
int grib_output_switch;
int netcdf_output_switch;
int pressure_level_output_switch;
int flight_level_output_switch;
int model_level_output_switch;
int surface_output_switch;
} Io_config;




