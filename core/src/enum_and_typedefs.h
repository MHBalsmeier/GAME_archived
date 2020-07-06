/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include <math.h>
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define M_V 0.0180152
#define R (N_A*K_B)
#define R_D (R/M_D)
#define R_V (R/M_V)
#define C_D_P 1005.0
#define C_D_V (C_D_P - R_D)
#define C_V_P 1858.0
#define C_V_V (C_V_P - R_V)
#define DELTA_C_V_P (C_V_P - C_D_P)
#define P_0 100000.0
#define RHO_WATER 1024.0
#define SECONDS_PER_MIN 60.0
#define SECONDS_PER_HOUR 3600.0
#define SECONDS_PER_DAY 86400.0
#define SECONDS_PER_YEAR 31536000.0
#define OMEGA (7.292115e-5)
#define SEMIMAJOR 6378137.0
#define SEMIMINOR 6356752.314
#define RADIUS pow(SEMIMAJOR*SEMIMAJOR*SEMIMINOR, 1.0/3.0)
#define H_BAR (1.054571817e-34)
#define entropy_constant_d (0.4*C_D_P*log(K_B/P_0*pow(M_D/N_A*K_B*exp(5.0/3)/(2*M_PI*H_BAR*H_BAR), 1.5)))
#define entropy_constant_v (0.4*C_V_P*log(K_B/P_0*pow(M_V/N_A*K_B*exp(5.0/3)/(2*M_PI*H_BAR*H_BAR), 1.5)))
#define EPSILON_TRACERS 0.00001

enum grid_integers {
RES_ID = 5,
NO_OF_TRACERS = 3,
NO_OF_CONDENSATED_TRACERS = 2,
NO_OF_SOLID_TRACERS = 1,
NO_OF_BASIC_TRIANGLES = 20,
NO_OF_PENTAGONS = 12,
NO_OF_HEXAGONS = (int) (10*(pow(2, 2*RES_ID) - 1)),
NO_OF_EDGES = 3*NO_OF_BASIC_TRIANGLES/2,
NO_OF_LAYERS = 26,
NO_OF_ORO_LAYERS = 17,
NO_OF_LEVELS = NO_OF_LAYERS + 1,
NO_OF_SCALARS_H = NO_OF_PENTAGONS + NO_OF_HEXAGONS,
NO_OF_VECTORS_H = (5*NO_OF_PENTAGONS/2 + 6/2*NO_OF_HEXAGONS),
NO_OF_VECTORS_V = NO_OF_SCALARS_H,
NO_OF_H_VECTORS = NO_OF_LAYERS*NO_OF_VECTORS_H,
NO_OF_V_VECTORS = NO_OF_LEVELS*NO_OF_VECTORS_V,
NO_OF_VECTORS_PER_LAYER = NO_OF_VECTORS_H + NO_OF_VECTORS_V,
NO_OF_TRIANGLES = (int) (NO_OF_BASIC_TRIANGLES*(pow(4, RES_ID))),
NO_OF_SCALARS = NO_OF_SCALARS_H*NO_OF_LAYERS,
NO_OF_VECTORS = NO_OF_H_VECTORS + NO_OF_V_VECTORS,
NO_OF_DUAL_SCALARS_H = NO_OF_TRIANGLES,
NO_OF_DUAL_VECTORS_H = 3*NO_OF_TRIANGLES/2,
NO_OF_DUAL_VECTORS_V = NO_OF_DUAL_SCALARS_H,
NO_OF_DUAL_H_VECTORS = NO_OF_LEVELS*NO_OF_DUAL_VECTORS_H,
NO_OF_DUAL_V_VECTORS = NO_OF_LAYERS*NO_OF_DUAL_VECTORS_V,
NO_OF_DUAL_VECTORS_PER_LAYER = NO_OF_DUAL_VECTORS_H + NO_OF_DUAL_VECTORS_V,
NO_OF_DUAL_SCALARS = NO_OF_LEVELS*NO_OF_DUAL_SCALARS_H,
NO_OF_DUAL_VECTORS = NO_OF_DUAL_H_VECTORS + NO_OF_DUAL_V_VECTORS,
TRIANGLES_PER_FACE = NO_OF_TRIANGLES/NO_OF_BASIC_TRIANGLES,
POINTS_PER_EDGE = (int) (pow(2, RES_ID) - 1),
SCALAR_POINTS_PER_INNER_FACE = (int) (0.5*(pow(2, RES_ID) - 2)*(pow(2, RES_ID) - 1)),
VECTOR_POINTS_PER_INNER_FACE = (int) (1.5*(pow(2, RES_ID) - 1)*pow(2, RES_ID))};

typedef double Scalar_field[NO_OF_SCALARS];
typedef double Vector_field[NO_OF_VECTORS];
typedef double Dual_vector_field[NO_OF_DUAL_VECTORS];
typedef double Curl_field[NO_OF_LAYERS*(NO_OF_DUAL_VECTORS_H + NO_OF_VECTORS_H) + NO_OF_DUAL_VECTORS_H];
typedef double Curl_field_one_layer[NO_OF_DUAL_VECTORS_H + NO_OF_VECTORS_H];
typedef double Tracer_densities[NO_OF_TRACERS*NO_OF_SCALARS];
typedef double Tracer_density_temperatures[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS];

typedef struct grid {
Vector_field normal_distance;
Scalar_field volume;
Vector_field area;
Scalar_field z_scalar;
Vector_field z_vector;
Scalar_field gravity_potential;
Vector_field gravity;
double e_kin_weights[12*NO_OF_SCALARS];
int e_kin_indices[12*NO_OF_SCALARS];
double vertical_contravar_unit[3*NO_OF_VECTORS_V*(NO_OF_ORO_LAYERS + 1)];
int from_index[NO_OF_VECTORS_H];
int to_index[NO_OF_VECTORS_H];
double direction[NO_OF_VECTORS_H];
int adjacent_vector_indices_h[6*NO_OF_SCALARS_H];
int adjacent_signs_h[6*NO_OF_SCALARS_H];
int recov_hor_par_curl_index[2*NO_OF_VECTORS_H];
double recov_hor_par_curl_weight[2*NO_OF_VECTORS_H];
int trsk_modified_velocity_indices[10*NO_OF_VECTORS_H];
int trsk_modified_curl_indices[10*NO_OF_VECTORS_H];
double trsk_modified_weights[10*NO_OF_VECTORS_H];
int recov_hor_ver_pri_index[4*NO_OF_VECTORS_H];
int recov_ver_index[6*NO_OF_VECTORS_V];
double recov_ver_0_pri_weight[6*NO_OF_VECTORS_V];
double recov_ver_0_curl_weight[6*NO_OF_VECTORS_V];
double recov_ver_1_pri_weight[6*NO_OF_VECTORS_V];
double recov_ver_1_curl_weight[6*NO_OF_VECTORS_V];
} Grid;

typedef struct dualgrid {
Curl_field_one_layer f_vec;
Curl_field area;
int vorticity_indices[4*NO_OF_VECTORS_H];
int vorticity_signs[4*NO_OF_VECTORS_H];
int h_curl_indices[4*NO_OF_DUAL_VECTORS_H];
int h_curl_signs[4*NO_OF_DUAL_VECTORS_H];
} Dualgrid;

typedef struct state {
Scalar_field entropy_gas;
Scalar_field density_dry;
Vector_field velocity_gas;
// density order: solid, liquid, vapour
Tracer_densities tracer_densities;
Tracer_density_temperatures tracer_density_temperatures;
} State;
