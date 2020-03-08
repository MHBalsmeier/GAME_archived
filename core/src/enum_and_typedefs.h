#include <math.h>
#define N_A 6.0221409e23
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define C_P 1005
#define C_V (C_P - R_D)
#define KAPPA (C_P/C_V)
#define P_0 100000
#define SECONDS_PER_HOUR 3600
#define SECONDS_PER_DAY 86400
#define SECONDS_PER_YEAR 31536000
#define OMEGA (7.292115e-5)
#define OMEGA_S (1.99099e-7)
#define SEMIMAJOR 6378137.0
#define SEMIMINOR 6356752.314
#define MASS (5.9723e24)
#define BETA 0.513812348
#define RHO_WAHTER 1024
#define SCALE_HEIGHT 8000

enum grid_integers {
RES_ID = 2,
NUMBER_OF_BASIC_TRIANGLES = 20,
NUMBER_OF_PENTAGONS = 12,
NUMBER_OF_HEXAGONS = (int) (10*(pow(2, 2*RES_ID) - 1)),
NUMBER_OF_EDGES = 3*NUMBER_OF_BASIC_TRIANGLES/2,
NUMBER_OF_LAYERS = (int) (3*pow(2, RES_ID - 2)),
NUMBER_OF_LEVELS = NUMBER_OF_LAYERS + 1,
NUMBER_OF_SCALARS_H = NUMBER_OF_PENTAGONS + NUMBER_OF_HEXAGONS,
NUMBER_OF_VECTORS_H = (5*NUMBER_OF_PENTAGONS/2 + 6/2*NUMBER_OF_HEXAGONS),
NUMBER_OF_H_VECTORS = NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H,
NUMBER_OF_V_VECTORS = NUMBER_OF_LEVELS*NUMBER_OF_SCALARS_H,
NUMBER_OF_TRIANGLES = (int) (NUMBER_OF_BASIC_TRIANGLES*(pow(4, RES_ID))),
NUMBER_OF_SCALARS = NUMBER_OF_SCALARS_H*NUMBER_OF_LAYERS,
NUMBER_OF_VECTORS = NUMBER_OF_H_VECTORS + NUMBER_OF_V_VECTORS,
NUMBER_OF_DUAL_SCALARS_H = NUMBER_OF_TRIANGLES,
NUMBER_OF_DUAL_VECTORS_H = 3*NUMBER_OF_TRIANGLES/2,
NUMBER_OF_DUAL_H_VECTORS = NUMBER_OF_LEVELS*NUMBER_OF_DUAL_VECTORS_H,
NUMBER_OF_DUAL_V_VECTORS = NUMBER_OF_LAYERS*NUMBER_OF_DUAL_SCALARS_H,
NUMBER_OF_DUAL_SCALARS = NUMBER_OF_LEVELS*NUMBER_OF_DUAL_SCALARS_H,
NUMBER_OF_DUAL_VECTORS = NUMBER_OF_DUAL_H_VECTORS + NUMBER_OF_DUAL_V_VECTORS,
TRIANGLES_PER_FACE = NUMBER_OF_TRIANGLES/NUMBER_OF_BASIC_TRIANGLES,
POINTS_PER_EDGE = (int) (pow(2, RES_ID) - 1),
SCALAR_POINTS_PER_INNER_FACE = (int) (0.5*(pow(2, RES_ID) - 2)*(pow(2, RES_ID) - 1)),
VECTOR_POINTS_PER_INNER_FACE = (int) (1.5*(pow(2, RES_ID) - 1)*pow(2, RES_ID))};

typedef double Scalar_field[NUMBER_OF_SCALARS];
typedef double Vector_field[NUMBER_OF_VECTORS];
typedef double Dual_scalar_field[NUMBER_OF_DUAL_SCALARS];
typedef double Dual_vector_field[NUMBER_OF_DUAL_VECTORS];

typedef struct grid {
Scalar_field latitude_scalar;
Scalar_field longitude_scalar;
Scalar_field z_scalar;
Vector_field latitude_vector;
Vector_field longitude_vector;
Vector_field z_vector;
Vector_field normal_distance;
Scalar_field volume;
Vector_field area;
double parallel_distance[2*NUMBER_OF_VECTORS];
double gravity[NUMBER_OF_V_VECTORS];
long from_indices[NUMBER_OF_VECTORS];
long to_indices[NUMBER_OF_VECTORS];
long adjacent_scalar_indices_h[6*NUMBER_OF_SCALARS];
long adjacent_scalar_index_lower[NUMBER_OF_SCALARS];
long adjacent_scalar_index_upper[NUMBER_OF_SCALARS];
long adjacent_vector_indices_h[6*NUMBER_OF_SCALARS];
long adjacent_vector_index_lower[NUMBER_OF_SCALARS];
long adjacent_vector_index_upper[NUMBER_OF_SCALARS];
short adjacent_signs_h[6*NUMBER_OF_SCALARS];
long vorticity_indices[6*NUMBER_OF_SCALARS];
short vorticity_signs[6*NUMBER_OF_SCALARS];
long h_curl_indices[4*NUMBER_OF_VECTORS];
short h_curl_signs[4*NUMBER_OF_VECTORS];
short vector_product_sign[NUMBER_OF_H_VECTORS];
long recov_hor_par_dual_index[11*NUMBER_OF_H_VECTORS];
double recov_hor_par_dual_weight[11*NUMBER_OF_H_VECTORS];
long recov_hor_ver_dual_index[2*NUMBER_OF_H_VECTORS];
double recov_hor_ver_dual_weight[2*NUMBER_OF_H_VECTORS];
long recov_hor_par_pri_index[2*NUMBER_OF_H_VECTORS];
double recov_hor_par_pri_weight[2*NUMBER_OF_H_VECTORS];
long recov_hor_ver_pri_index[4*NUMBER_OF_H_VECTORS];
double recov_hor_ver_pri_weight[4*NUMBER_OF_H_VECTORS];
long recov_ver_1_pri_index[6*NUMBER_OF_V_VECTORS];
double recov_ver_1_pri_weight[6*NUMBER_OF_V_VECTORS];
long recov_ver_1_dual_index[6*NUMBER_OF_V_VECTORS];
double recov_ver_1_dual_weight[6*NUMBER_OF_V_VECTORS];
long recov_ver_2_pri_index[6*NUMBER_OF_V_VECTORS];
double recov_ver_2_pri_weight[6*NUMBER_OF_V_VECTORS];
long recov_ver_2_dual_index[6*NUMBER_OF_V_VECTORS];
double recov_ver_2_dual_weight[6*NUMBER_OF_V_VECTORS];
} Grid;

typedef struct dualgrid {
Dual_scalar_field latitude_scalar;
Dual_scalar_field longitude_scalar;
Dual_scalar_field z_scalar;
Dual_vector_field latitude_vector;
Dual_vector_field longitude_vector;
Dual_vector_field z_vector;
Dual_vector_field normal_distance;
Dual_vector_field area;
Dual_vector_field f_vec;
double parallel_distance[2*NUMBER_OF_DUAL_VECTORS];
long to_indices[NUMBER_OF_DUAL_VECTORS];
long from_indices[NUMBER_OF_DUAL_VECTORS];
long adjacent_scalar_indices_h[3*NUMBER_OF_DUAL_SCALARS];
long adjacent_scalar_index_lower[NUMBER_OF_DUAL_SCALARS];
long adjacent_scalar_index_upper[NUMBER_OF_DUAL_SCALARS];
long adjacent_vector_indices_h[3*NUMBER_OF_DUAL_SCALARS];
long adjacent_vector_index_lower[NUMBER_OF_DUAL_SCALARS];
long adjacent_vector_index_upper[NUMBER_OF_DUAL_SCALARS];
short adjacent_signs_h[3*NUMBER_OF_DUAL_SCALARS];
long vorticity_indices[3*NUMBER_OF_DUAL_V_VECTORS];
short vorticity_signs[3*NUMBER_OF_DUAL_V_VECTORS];
long h_curl_indices[4*NUMBER_OF_DUAL_H_VECTORS];
short h_curl_signs[4*NUMBER_OF_DUAL_H_VECTORS];
} Dualgrid;

typedef struct state {
Scalar_field density;
Scalar_field pot_temp;
Vector_field wind;
} State;
