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
#define SEMIMAJOR 6378137.0
#define SEMIMINOR 6356752.314

enum grid_integers {
RES_ID = 5,
NUMBER_OF_BASIC_TRIANGLES = 20,
NUMBER_OF_PENTAGONS = 12,
NUMBER_OF_HEXAGONS = (int) (10*(pow(2, 2*RES_ID) - 1)),
NUMBER_OF_EDGES = 3*NUMBER_OF_BASIC_TRIANGLES/2,
NUMBER_OF_LAYERS = 26,
NUMBER_OF_ORO_LAYERS = 17,
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