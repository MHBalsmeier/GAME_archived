/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "geos95.h"

int tangential_wind(Vector_field in_field, int layer_index, int h_index, double *component, Grid *grid)
{
    *component = 0;
	for (int i = 0; i < 10; ++i)
	{
		*component += grid -> trsk_weights[10*h_index + i]*in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> trsk_indices[10*h_index + i]];
	}
    return 0;
}

int calc_uv_at_edge(Vector_field in_field, Vector_field out_field_u, Vector_field out_field_v, Grid *grid)
{
	/*
	This function diagnozes eastward and northward components of a vector field at edges.
	*/
	int layer_index, h_index;
	double wind_0, wind_1;
	#pragma omp parallel for private(layer_index, h_index, wind_0, wind_1)
	for (int i = 0; i < NO_OF_H_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		wind_0 = in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
		tangential_wind(in_field, layer_index, h_index, &wind_1, grid);
		passive_turn(wind_0, wind_1, -grid -> direction[h_index],
		&out_field_u[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index],
		&out_field_v[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index]);
    }
	return 0;
}

