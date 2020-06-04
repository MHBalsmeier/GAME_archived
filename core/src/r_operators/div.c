#include "../enum_and_typedefs.h"
#include "../diagnostics/diagnostics.h"
#include <stdio.h>

int divergence(Vector_field in_field, Scalar_field out_field, Grid *grid, int allow_surface_flux)
{
    int number_of_edges, layer_index, h_index, retval;
    double contra_upper, contra_lower, comp_h, comp_v;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        comp_h = 0;
        layer_index = i/NUMBER_OF_SCALARS_H;
        h_index = i - layer_index*NUMBER_OF_SCALARS_H;
        number_of_edges = 6;
        if (h_index < NUMBER_OF_PENTAGONS)
            number_of_edges = 5;
        for (int j = 0; j < number_of_edges; ++j)
            comp_h += in_field[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]*grid -> adjacent_signs_h[6*h_index + j]*grid -> area[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]];
        if (layer_index == 0)
        {
        	retval = vertical_contravariant_normalized(in_field, layer_index + 1, h_index, grid, &contra_lower);
            if (retval != 0)
            	printf("Error in vertical_contravariant_normalized called at position 0 from divergence.\n");
            comp_v = -contra_lower*grid -> area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        }
        else if (layer_index == NUMBER_OF_LAYERS - 1)
        {
        	retval = vertical_contravariant_normalized(in_field, layer_index, h_index, grid, &contra_upper);
            if (retval != 0)
            	printf("Error in vertical_contravariant_normalized called at position 1 from divergence.\n");
        	retval = vertical_contravariant_normalized(in_field, layer_index + 1, h_index, grid, &contra_lower);
            if (retval != 0)
            	printf("Error in vertical_contravariant_normalized called at position 2 from divergence.\n");
            comp_v = contra_upper*grid -> area[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER] - allow_surface_flux*contra_lower*grid -> area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        }
        else
        {
        	retval = vertical_contravariant_normalized(in_field, layer_index, h_index, grid, &contra_upper);
            if (retval != 0)
            	printf("Error in vertical_contravariant_normalized called at position 3 from divergence.\n");
        	retval = vertical_contravariant_normalized(in_field, layer_index + 1, h_index, grid, &contra_lower);
            if (retval != 0)
            	printf("Error in vertical_contravariant_normalized called at position 4 from divergence.\n");
            comp_v = contra_upper*grid -> area[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER] - contra_lower*grid -> area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        }
        out_field[i] = 1/grid -> volume[i]*(comp_h + comp_v);
    }
    return 0;
}
