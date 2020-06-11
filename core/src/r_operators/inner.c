#include "../enum_and_typedefs.h"
#include <stdio.h>

int inner(Vector_field in_field_0, Vector_field in_field_1, Scalar_field out_field, Grid *grid, Dualgrid *dualgrid)
{
    int layer_index, h_index, number_of_edges, vector_index, dual_vector_index;
    double comp_h, comp_v, area_pre_factor, normal_distance_dual;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
        h_index = i - layer_index*NUMBER_OF_SCALARS_H;
		number_of_edges = 6;
		if (h_index < NUMBER_OF_PENTAGONS)
			number_of_edges = 5;		
		comp_h = 0;
		area_pre_factor = 0.25/(grid -> area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER]*pow((RADIUS + grid -> z_scalar[i])/(RADIUS + grid -> z_vector[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER]), 2));
		for (int j = 0; j < number_of_edges; ++j)
		{
			vector_index = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j];
			dual_vector_index = grid -> adjacent_vector_indices_h[6*h_index + j] + (layer_index + 1)*NUMBER_OF_DUAL_VECTORS_PER_LAYER;
			normal_distance_dual = dualgrid -> normal_distance[dual_vector_index]*(RADIUS + grid -> z_scalar[i])/(RADIUS + grid -> z_vector[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER]);
			comp_h += area_pre_factor*grid -> normal_distance[vector_index]*normal_distance_dual*in_field_0[vector_index]*in_field_1[vector_index];
		}
        comp_v = in_field_0[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER]*in_field_1[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER];
        comp_v += in_field_0[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER]*in_field_1[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        out_field[i] = comp_h + 0.5*comp_v;
    }
    return 0;
}
