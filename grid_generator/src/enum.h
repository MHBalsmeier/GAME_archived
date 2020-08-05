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
RES_ID = 6,
NO_OF_BASIC_TRIANGLES = 20,
NO_OF_PENTAGONS = 12,
NO_OF_HEXAGONS = (int) (10*(pow(2, 2*RES_ID) - 1)),
NO_OF_EDGES = 3*NO_OF_BASIC_TRIANGLES/2,
NO_OF_LAYERS = 26,
NO_OF_ORO_LAYERS = 17,
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
