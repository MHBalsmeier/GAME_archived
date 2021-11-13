/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, constants and types are defined.
*/

#include <math.h>

enum grid_integers {
// This determines the horizontal resolution.
RES_ID = 5,
// This has to conform with the grid file and the initialization state file.
NO_OF_LAYERS = 26,
// The number of layers affected by orography. This also has to conform with the grid file and the initialization state file.
NO_OF_GASEOUS_CONSTITUENTS = 2,
NO_OF_CONDENSED_CONSTITUENTS = 4,
// the number of blocks into which the arrays will be split up for the radiation calculation
// (NO_OF_SCALARS_H must be divisible by this number)
NO_OF_RAD_BLOCKS = 18,
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
NO_OF_SCALARS_RAD = NO_OF_SCALARS/NO_OF_RAD_BLOCKS,
NO_OF_SCALARS_RAD_PER_LAYER = NO_OF_SCALARS_RAD/NO_OF_LAYERS,
NO_OF_DUAL_SCALARS_H = NO_OF_TRIANGLES,
NO_OF_DUAL_H_VECTORS = NO_OF_LEVELS*NO_OF_VECTORS_H,
NO_OF_DUAL_V_VECTORS = NO_OF_LAYERS*NO_OF_DUAL_SCALARS_H,
NO_OF_DUAL_VECTORS_PER_LAYER = NO_OF_VECTORS_H + NO_OF_DUAL_SCALARS_H,
NO_OF_DUAL_SCALARS = NO_OF_LEVELS*NO_OF_DUAL_SCALARS_H,
NO_OF_DUAL_VECTORS = NO_OF_DUAL_H_VECTORS + NO_OF_DUAL_V_VECTORS,
NO_OF_LON_IO_POINTS = (int) (4*(pow(2, RES_ID))),
NO_OF_LAT_IO_POINTS = (int) (2*(pow(2, RES_ID))),
NO_OF_LATLON_IO_POINTS = NO_OF_LON_IO_POINTS*NO_OF_LAT_IO_POINTS,
POINTS_PER_EDGE = (int) (pow(2, RES_ID) - 1),
TRIANGLES_PER_FACE = NO_OF_TRIANGLES/NO_OF_BASIC_TRIANGLES,
SCALAR_POINTS_PER_INNER_FACE = (int) (0.5*(pow(2, RES_ID) - 2)*(pow(2, RES_ID) - 1)),
VECTOR_POINTS_PER_INNER_FACE = (int) (1.5*(pow(2, RES_ID) - 1)*pow(2, RES_ID))};

typedef double Scalar_field[NO_OF_SCALARS];
typedef double Vector_field[NO_OF_VECTORS];
typedef double Dual_vector_field[NO_OF_DUAL_VECTORS];
typedef double Curl_field[NO_OF_LAYERS*2*NO_OF_VECTORS_H + NO_OF_VECTORS_H];
// all constituents have a mass density
typedef double Mass_densities[NO_OF_CONSTITUENTS*NO_OF_SCALARS];
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
Scalar_field theta_bg;
Scalar_field exner_bg;
Vector_field exner_bg_grad;
int trsk_indices[10*NO_OF_VECTORS_H];
int trsk_modified_curl_indices[10*NO_OF_VECTORS_H];
int from_index[NO_OF_VECTORS_H];
int to_index[NO_OF_VECTORS_H];
int adjacent_vector_indices_h[6*NO_OF_SCALARS_H];
int adjacent_signs_h[6*NO_OF_SCALARS_H];
int density_to_rhombi_indices[4*NO_OF_VECTORS_H];
int no_of_shaded_points_scalar[NO_OF_SCALARS_H];
int no_of_shaded_points_vector[NO_OF_VECTORS_H];
double latitude_scalar[NO_OF_SCALARS_H];
double longitude_scalar[NO_OF_SCALARS_H];
double inner_product_weights[8*NO_OF_SCALARS];
double direction[NO_OF_VECTORS_H];
double density_to_rhombi_weights[4*NO_OF_VECTORS_H];
double trsk_weights[10*NO_OF_VECTORS_H];
double sfc_albedo[NO_OF_SCALARS_H];
double sfc_rho_c[NO_OF_SCALARS_H];
int is_land[NO_OF_SCALARS_H];
int latlon_interpol_indices[5*NO_OF_LATLON_IO_POINTS];
double latlon_interpol_weights[5*NO_OF_LATLON_IO_POINTS];
double stretching_parameter;
double mean_area_cell;
} Grid;

// Contains properties of the dual grid.
typedef struct dualgrid {
Curl_field area;
Dual_vector_field z_vector;
Dual_vector_field normal_distance;
int from_index[NO_OF_VECTORS_H];
int to_index[NO_OF_VECTORS_H];
int vorticity_indices_triangles[3*NO_OF_DUAL_SCALARS_H];
int vorticity_signs_triangles[3*NO_OF_DUAL_SCALARS_H];
double f_vec[2*NO_OF_VECTORS_H];
} Dualgrid;

typedef struct state {
Mass_densities rho; // density order: solid, liquid, vapour
Scalar_field rhotheta;
Scalar_field theta_pert;
Scalar_field exner_pert;
Condensed_density_temperatures condensed_density_temperatures;
Vector_field wind;
} State;

// Collects diagnostic quantities. Note: in fact, forcings are also diagnostic quantities.
typedef struct diagnostics {
Vector_field flux_density;
Scalar_field flux_density_divv;
double rel_vort_on_triangles[NO_OF_DUAL_V_VECTORS];
Curl_field rel_vort;
Curl_field pot_vort;
Scalar_field temperature_gas;
Scalar_field c_g_p_field;
Scalar_field v_squared;
Scalar_field wind_divv;
Vector_field curl_of_vorticity;
Scalar_field scalar_field_placeholder;
Vector_field vector_field_placeholder;
Vector_field u_at_edge;
Vector_field v_at_edge;
Scalar_field u_at_cell;
Scalar_field v_at_cell;
double dv_hdz[NO_OF_H_VECTORS + NO_OF_VECTORS_H];
} Diagnostics;

// needed for the radiation calculation
typedef struct radiation {
double lat_scal[NO_OF_SCALARS_RAD_PER_LAYER];
double lon_scal[NO_OF_SCALARS_RAD_PER_LAYER];
double sfc_sw_in[NO_OF_SCALARS_RAD_PER_LAYER];
double sfc_lw_out[NO_OF_SCALARS_RAD_PER_LAYER];
double sfc_albedo[NO_OF_SCALARS_RAD_PER_LAYER];
double temp_sfc[NO_OF_SCALARS_RAD_PER_LAYER];
double z_scal[NO_OF_SCALARS_RAD];
double z_vect[NO_OF_SCALARS_RAD + NO_OF_SCALARS_RAD_PER_LAYER];
double rho[NO_OF_CONSTITUENTS*NO_OF_SCALARS_RAD];
double temp[NO_OF_SCALARS_RAD];
double rad_tend[NO_OF_SCALARS_RAD];
} Radiation;

// Collects forcings.
typedef struct forcings {
Vector_field pgrad_acc_old;
Vector_field pressure_gradient_acc_neg_nl;
Vector_field pressure_gradient_acc_neg_l;
Vector_field v_squared_grad;
Vector_field pot_vort_tend;
double sfc_sw_in[NO_OF_SCALARS_H];
double sfc_lw_out[NO_OF_SCALARS_H];
Scalar_field radiation_tendency;
} Forcings;

// Info on the run configuration is collected here.
typedef struct config_info {
int totally_first_step_bool;
int temperature_diff_h;
int temperature_diff_v;
int momentum_diff_h;
int momentum_diff_v;
int tracer_diff_h;
int tracer_diff_v;
int rad_on;
int rad_update;
int assume_lte;
int slow_fast_ratio;
int nwp_mode;
int delta_t_between_analyses;
int explicit_boundary_layer;
int waves_on;
double diff_h_smag_fac;
double shear_bg;
double damping_start_height_over_toa;
double damping_coeff_max;
} Config_info;

// soil information
typedef struct soil {
double temperature[NO_OF_SCALARS_H];
double power_flux_density_sensible[NO_OF_SCALARS_H];
} Soil;

// Contains everything on turbulence parametrizations as well as constituent-related quantities.
typedef struct irreversible_quantities {
Scalar_field temperature_diffusion_heating;
Vector_field friction_acc;
Scalar_field heating_diss;
Scalar_field scalar_diffusion_coeff_numerical_h;
Scalar_field scalar_diffusion_coeff_numerical_v;
Scalar_field pressure_gradient_decel_factor;
double mass_source_rates[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS];
double constituent_heat_source_rates[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS];
Scalar_field viscosity_div_eff;
Vector_field viscosity_curl_eff_rhombi;
double viscosity_curl_eff_triangles[NO_OF_DUAL_V_VECTORS];
double vert_hor_viscosity_eff[NO_OF_H_VECTORS + NO_OF_VECTORS_H];
Scalar_field tke;
} Irreversible_quantities;

// Info on input and output is collected here.
typedef struct io_config {
int grib_output_switch;
int netcdf_output_switch;
int pressure_level_output_switch;
int model_level_output_switch;
int surface_output_switch;
} Io_config;




