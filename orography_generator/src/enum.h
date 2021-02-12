/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include <math.h>

enum grid_integers {
RES_ID = 5,
NO_OF_BASIC_TRIANGLES = 20,
NO_OF_PENTAGONS = 12,
NO_OF_HEXAGONS = (int) (10*(pow(2, 2*RES_ID) - 1)),
NO_OF_EDGES = 3*NO_OF_BASIC_TRIANGLES/2,
NO_OF_LAYERS = 6,
NO_OF_ORO_LAYERS = 4,
NO_OF_LEVELS = NO_OF_LAYERS + 1,
NO_OF_SCALARS_H = NO_OF_PENTAGONS + NO_OF_HEXAGONS,
NO_OF_VECTORS_H = (5*NO_OF_PENTAGONS/2 + 6/2*NO_OF_HEXAGONS),
NO_OF_VECTORS_PER_LAYER = NO_OF_VECTORS_H + NO_OF_SCALARS_H,
NO_OF_TRIANGLES = (int) (NO_OF_BASIC_TRIANGLES*(pow(4, RES_ID))),
NO_OF_DUAL_SCALARS_H = NO_OF_TRIANGLES,
TRIANGLES_PER_FACE = NO_OF_TRIANGLES/NO_OF_BASIC_TRIANGLES,
POINTS_PER_EDGE = (int) (pow(2, RES_ID) - 1),
SCALAR_POINTS_PER_INNER_FACE = (int) (0.5*(pow(2, RES_ID) - 2)*(pow(2, RES_ID) - 1)),
VECTOR_POINTS_PER_INNER_FACE = (int) (1.5*(pow(2, RES_ID) - 1)*pow(2, RES_ID))};
