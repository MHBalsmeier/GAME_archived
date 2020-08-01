/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "../diagnostics/diagnostics.h"
#include <stdio.h>

int divv_h(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
    int number_of_edges, layer_index, h_index;
    double contra_upper, contra_lower, comp_h, comp_v;
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        comp_h = 0;
        layer_index = i/NO_OF_SCALARS_H;
        h_index = i - layer_index*NO_OF_SCALARS_H;
        number_of_edges = 6;
        if (h_index < NO_OF_PENTAGONS)
        {
            number_of_edges = 5;
        }
        for (int j = 0; j < number_of_edges; ++j)
        {
            comp_h += in_field[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]*grid -> adjacent_signs_h[6*h_index + j]*grid -> area[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]];
        }
        if (layer_index == 0)
        {
            comp_v = 0;
        }
        else if (layer_index == NO_OF_LAYERS - 1)
        {
        	vertical_contravariant_normalized_h(in_field, layer_index, h_index, grid, &contra_upper);
            comp_v = contra_upper*grid -> area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
        }
        else
        {
            contra_upper = 0;
            if (layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
            {
                vertical_contravariant_normalized_h(in_field, layer_index, h_index, grid, &contra_upper);
            }
            contra_lower = 0;
            if (layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS - 1)
            {
                vertical_contravariant_normalized_h(in_field, layer_index + 1, h_index, grid, &contra_lower);
            }
            comp_v = contra_upper*grid -> area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER] - contra_lower*grid -> area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        }
        out_field[i] = 1/grid -> volume[i]*(comp_h + comp_v);
    }
    return 0;
}

int divv_v_columns(double in_vector[], double out_vector[], int h_index, Grid *grid)
{
	// In_vector must be contravariant already.
    double contra_upper, contra_lower, comp_v;
    for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
    {
        if (layer_index == 0)
        {
            contra_lower = in_vector[0];
            comp_v = -contra_lower*grid -> area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        }
        else if (layer_index == NO_OF_LAYERS - 1)
        {
        	contra_upper = in_vector[layer_index - 1];
            comp_v = contra_upper*grid -> area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
        }
        else
        {
            contra_upper = in_vector[layer_index - 1];
            contra_lower = in_vector[layer_index];
            comp_v = contra_upper*grid -> area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER] - contra_lower*grid -> area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        }
        out_vector[layer_index] = 1/grid -> volume[h_index + layer_index*NO_OF_SCALARS_H]*comp_v;
    }
    return 0;
}







