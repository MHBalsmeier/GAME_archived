#include <math.h>
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define C_P 1005.0
#define C_V (C_P - R_D)
#define KAPPA (C_P/C_V)
#define P_0 100000.0
#define SECONDS_PER_MIN 60.0
#define SECONDS_PER_HOUR 3600.0
#define SECONDS_PER_DAY 86400.0
#define SECONDS_PER_YEAR 31536000.0
#define OMEGA (7.292115e-5)
#define OMEGA_S (1.99099e-7)
#define SEMIMAJOR 6378137.0
#define SEMIMINOR 6356752.314
#define MASS (5.9723e24)
#define BETA 0.513812348
#define RHO_WAHTER 1024.0
#define SCALE_HEIGHT 8000.0

enum grid_integers {
RES_ID = 4,
NUMBER_OF_BASIC_TRIANGLES = 20,
NUMBER_OF_PENTAGONS = 12,
NUMBER_OF_HEXAGONS = (int) (10*(pow(2, 2*RES_ID) - 1)),
NUMBER_OF_EDGES = 3*NUMBER_OF_BASIC_TRIANGLES/2,
NUMBER_OF_LAYERS = (int) fmax((3*pow(2, RES_ID - 3)), 3),
NUMBER_OF_LEVELS = NUMBER_OF_LAYERS + 1,
NUMBER_OF_SCALARS_H = NUMBER_OF_PENTAGONS + NUMBER_OF_HEXAGONS,
NUMBER_OF_VECTORS_H = (5*NUMBER_OF_PENTAGONS/2 + 6/2*NUMBER_OF_HEXAGONS),
NUMBER_OF_VECTORS_V = NUMBER_OF_SCALARS_H,
NUMBER_OF_H_VECTORS = NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H,
NUMBER_OF_V_VECTORS = NUMBER_OF_LEVELS*NUMBER_OF_SCALARS_H,
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

typedef struct grid {
Vector_field gravity;
Vector_field normal_distance;
Scalar_field volume;
Vector_field area;
long from_index[NUMBER_OF_VECTORS_H];
long to_index[NUMBER_OF_VECTORS_H];
long adjacent_vector_indices_h[6*NUMBER_OF_SCALARS_H];
short adjacent_signs_h[6*NUMBER_OF_SCALARS_H];
long vorticity_indices[6*NUMBER_OF_SCALARS_H];
short vorticity_signs[6*NUMBER_OF_SCALARS_H];
long h_curl_indices[4*NUMBER_OF_VECTORS_H];
short h_curl_signs[4*NUMBER_OF_VECTORS_H];
short vector_product_sign[NUMBER_OF_VECTORS_H];
long recov_hor_par_dual_index[2*NUMBER_OF_VECTORS_H];
double recov_hor_par_dual_weight[2*NUMBER_OF_VECTORS_H];
long recov_hor_ver_dual_index[2*NUMBER_OF_VECTORS_H];
double recov_hor_ver_dual_weight[2*NUMBER_OF_VECTORS_H];
long recov_hor_par_pri_index[2*NUMBER_OF_VECTORS_H];
double recov_hor_par_pri_weight[2*NUMBER_OF_VECTORS_H];
long recov_hor_ver_pri_index[4*NUMBER_OF_VECTORS_H];
double recov_hor_ver_pri_weight[4*NUMBER_OF_VECTORS_H];
long recov_ver_0_pri_index[6*NUMBER_OF_VECTORS_V];
double recov_ver_0_pri_weight[6*NUMBER_OF_VECTORS_V];
long recov_ver_0_dual_index[6*NUMBER_OF_VECTORS_V];
double recov_ver_0_dual_weight[6*NUMBER_OF_VECTORS_V];
long recov_ver_1_pri_index[6*NUMBER_OF_VECTORS_V];
double recov_ver_1_pri_weight[6*NUMBER_OF_VECTORS_V];
long recov_ver_1_dual_index[6*NUMBER_OF_VECTORS_V];
double recov_ver_1_dual_weight[6*NUMBER_OF_VECTORS_V];
} Grid;

typedef struct dualgrid {
Dual_vector_field normal_distance;
Dual_vector_field area;
double f_vec[NUMBER_OF_DUAL_VECTORS_PER_LAYER];
long to_index[NUMBER_OF_DUAL_VECTORS_H];
long from_index[NUMBER_OF_DUAL_VECTORS_H];
long vorticity_indices[3*NUMBER_OF_DUAL_VECTORS_V];
short vorticity_signs[3*NUMBER_OF_DUAL_VECTORS_V];
long h_curl_indices[4*NUMBER_OF_DUAL_VECTORS_H];
short h_curl_signs[4*NUMBER_OF_DUAL_VECTORS_H];
} Dualgrid;

typedef struct state {
Scalar_field density_pot_temp;
Scalar_field density;
Vector_field wind;
} State;
