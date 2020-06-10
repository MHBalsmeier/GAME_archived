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
RES_ID = 4,
NUMBER_OF_ADD_COMPS = 3,
NUMBER_OF_COND_ADD_COMPS = 2,
NUMBER_OF_SOLID_ADD_COMPS = 1,
NUMBER_OF_BASIC_TRIANGLES = 20,
NUMBER_OF_PENTAGONS = 12,
NUMBER_OF_HEXAGONS = (int) (10*(pow(2, 2*RES_ID) - 1)),
NUMBER_OF_EDGES = 3*NUMBER_OF_BASIC_TRIANGLES/2,
NUMBER_OF_LAYERS = 6,
NUMBER_OF_ORO_LAYERS = 4,
NUMBER_OF_LEVELS = NUMBER_OF_LAYERS + 1,
NUMBER_OF_SCALARS_H = NUMBER_OF_PENTAGONS + NUMBER_OF_HEXAGONS,
NUMBER_OF_VECTORS_H = (5*NUMBER_OF_PENTAGONS/2 + 6/2*NUMBER_OF_HEXAGONS),
NUMBER_OF_VECTORS_V = NUMBER_OF_SCALARS_H,
NUMBER_OF_H_VECTORS = NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H,
NUMBER_OF_V_VECTORS = NUMBER_OF_LEVELS*NUMBER_OF_VECTORS_V,
NUMBER_OF_VECTORS_PER_LAYER = NUMBER_OF_VECTORS_H + NUMBER_OF_VECTORS_V,
NUMBER_OF_TRIANGLES = (int) (NUMBER_OF_BASIC_TRIANGLES*(pow(4, RES_ID))),
NUMBER_OF_SCALARS = NUMBER_OF_SCALARS_H*NUMBER_OF_LAYERS,
NUMBER_OF_VECTORS = NUMBER_OF_H_VECTORS + NUMBER_OF_V_VECTORS,
NUMBER_OF_DUAL_SCALARS_H = NUMBER_OF_TRIANGLES,
NUMBER_OF_DUAL_VECTORS_H = 3*NUMBER_OF_TRIANGLES/2,
NUMBER_OF_DUAL_VECTORS_V = NUMBER_OF_DUAL_SCALARS_H,
NUMBER_OF_DUAL_H_VECTORS = NUMBER_OF_LEVELS*NUMBER_OF_DUAL_VECTORS_H,
NUMBER_OF_DUAL_V_VECTORS = NUMBER_OF_LAYERS*NUMBER_OF_DUAL_VECTORS_V,
NUMBER_OF_DUAL_VECTORS_PER_LAYER = NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_DUAL_VECTORS_V,
NUMBER_OF_DUAL_SCALARS = NUMBER_OF_LEVELS*NUMBER_OF_DUAL_SCALARS_H,
NUMBER_OF_DUAL_VECTORS = NUMBER_OF_DUAL_H_VECTORS + NUMBER_OF_DUAL_V_VECTORS,
TRIANGLES_PER_FACE = NUMBER_OF_TRIANGLES/NUMBER_OF_BASIC_TRIANGLES,
POINTS_PER_EDGE = (int) (pow(2, RES_ID) - 1),
SCALAR_POINTS_PER_INNER_FACE = (int) (0.5*(pow(2, RES_ID) - 2)*(pow(2, RES_ID) - 1)),
VECTOR_POINTS_PER_INNER_FACE = (int) (1.5*(pow(2, RES_ID) - 1)*pow(2, RES_ID))};

typedef double Scalar_field[NUMBER_OF_SCALARS];
typedef double Vector_field[NUMBER_OF_VECTORS];
typedef double Dual_vector_field[NUMBER_OF_DUAL_VECTORS];
typedef double Add_comp_densities[NUMBER_OF_ADD_COMPS*NUMBER_OF_SCALARS];
typedef double Add_comp_temps[NUMBER_OF_COND_ADD_COMPS*NUMBER_OF_SCALARS];

typedef struct grid {
Vector_field normal_distance;
Scalar_field volume;
Vector_field area;
Scalar_field z_scalar;
Vector_field z_vector;
Scalar_field gravity_potential;
Vector_field gravity;
double vertical_contravar_unit[3*NUMBER_OF_VECTORS_V*(NUMBER_OF_ORO_LAYERS + 1)];
int from_index[NUMBER_OF_VECTORS_H];
int to_index[NUMBER_OF_VECTORS_H];
double direction[NUMBER_OF_VECTORS_H];
double z_surface[NUMBER_OF_SCALARS_H];
int adjacent_vector_indices_h[6*NUMBER_OF_SCALARS_H];
int adjacent_signs_h[6*NUMBER_OF_SCALARS_H];
int recov_hor_par_dual_index[2*NUMBER_OF_VECTORS_H];
double recov_hor_par_dual_weight[2*NUMBER_OF_VECTORS_H];
int recov_hor_ver_dual_index[2*NUMBER_OF_VECTORS_H];
double recov_hor_ver_dual_weight[2*NUMBER_OF_VECTORS_H];
int recov_hor_par_pri_index[10*NUMBER_OF_VECTORS_H];
double recov_hor_par_pri_weight[10*NUMBER_OF_VECTORS_H];
int recov_hor_ver_pri_index[4*NUMBER_OF_VECTORS_H];
double recov_hor_ver_pri_weight[4*NUMBER_OF_VECTORS_H];
int recov_ver_0_pri_index[6*NUMBER_OF_VECTORS_V];
double recov_ver_0_pri_weight[6*NUMBER_OF_VECTORS_V];
int recov_ver_0_dual_index[6*NUMBER_OF_VECTORS_V];
double recov_ver_0_dual_weight[6*NUMBER_OF_VECTORS_V];
int recov_ver_1_pri_index[6*NUMBER_OF_VECTORS_V];
double recov_ver_1_pri_weight[6*NUMBER_OF_VECTORS_V];
int recov_ver_1_dual_index[6*NUMBER_OF_VECTORS_V];
double recov_ver_1_dual_weight[6*NUMBER_OF_VECTORS_V];
} Grid;

typedef struct dualgrid {
Dual_vector_field normal_distance;
Dual_vector_field z_vector;
double area[NUMBER_OF_DUAL_H_VECTORS + NUMBER_OF_H_VECTORS];
double f_vec[NUMBER_OF_DUAL_VECTORS_PER_LAYER];
int vorticity_indices[4*NUMBER_OF_VECTORS_H];
int vorticity_signs[4*NUMBER_OF_VECTORS_H];
int h_curl_indices[4*NUMBER_OF_DUAL_VECTORS_H];
int h_curl_signs[4*NUMBER_OF_DUAL_VECTORS_H];
int adjacent_vector_indices_h[3*NUMBER_OF_DUAL_SCALARS_H];
} Dualgrid;

typedef struct state {
Scalar_field entropy;
Scalar_field density;
Vector_field wind;
// density order: solid, liquid, vapour
Add_comp_densities add_comp_densities;
Add_comp_temps add_comp_temps;
} State;
