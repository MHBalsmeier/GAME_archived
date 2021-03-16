/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>

int scalar_times_vector(Scalar_field scalar_field, Vector_field vector_field, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, lower_index, upper_index, i;
    double scalar_value;
	#pragma omp parallel for private (layer_index, h_index, lower_index, upper_index, i, scalar_value)
    for (i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            scalar_value
            = 0.5*(
            scalar_field[grid -> to_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]
            + scalar_field[grid -> from_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]);
        }
        else
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            scalar_value = 0.5*(
            scalar_field[upper_index]
            + scalar_field[lower_index]);
        }
    	out_field[i] = scalar_value*vector_field[i];
    }
    // linear extrapolation to the TOA
    #pragma omp parallel for private(scalar_value)
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        scalar_value
        = scalar_field[i]
        + (scalar_field[i] - scalar_field[i + NO_OF_SCALARS_H])
        /(grid -> z_scalar[i] - grid -> z_scalar[i + NO_OF_SCALARS_H])
        *(grid -> z_vector[i] - grid -> z_scalar[i]);
    	out_field[i] = scalar_value*vector_field[i];
    }
    // linear extrapolation to the surface
    #pragma omp parallel for private(layer_index, h_index, upper_index, scalar_value)
    for (int i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
        scalar_value
        = scalar_field[upper_index]
        + (scalar_field[upper_index - NO_OF_SCALARS_H] - scalar_field[upper_index])
        /(grid -> z_scalar[upper_index - NO_OF_SCALARS_H] - grid -> z_scalar[upper_index])
        *(grid -> z_vector[i] - grid -> z_scalar[upper_index]);
    	out_field[i] = scalar_value*vector_field[i];
    }
    return 0;
}

int scalar_times_vector_scalar_h_v(Scalar_field in_field_h, Scalar_field in_field_v, Vector_field vector_field, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, lower_index, upper_index, i;
    double scalar_value;
    #pragma omp parallel for private (layer_index, h_index, lower_index, upper_index, i, scalar_value)
    for (i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            scalar_value
            = 0.5*(
            in_field_h[grid -> to_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]
            + in_field_h[grid -> from_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]);
        }
        else
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            scalar_value = 0.5*(
            in_field_v[upper_index]
            + in_field_v[lower_index]);
        }
        out_field[i] = scalar_value*vector_field[i];
    }
    // linear extrapolation to the TOA
    #pragma omp parallel for private(scalar_value)
    for (i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        scalar_value
        = in_field_v[i]
        + (in_field_v[i] - in_field_v[i + NO_OF_SCALARS_H])
        /(grid -> z_scalar[i] - grid -> z_scalar[i + NO_OF_SCALARS_H])
        *(grid -> z_vector[i] - grid -> z_scalar[i]);
        out_field[i] = scalar_value*vector_field[i];
    }
    // linear extrapolation to the surface
    #pragma omp parallel for private(layer_index, h_index, upper_index, scalar_value)
    for (i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
        scalar_value
        = in_field_v[upper_index]
        + (in_field_v[upper_index - NO_OF_SCALARS_H] - in_field_v[upper_index])
        /(grid -> z_scalar[upper_index - NO_OF_SCALARS_H] - grid -> z_scalar[upper_index])
        *(grid -> z_vector[i] - grid -> z_scalar[upper_index]);
        out_field[i] = scalar_value*vector_field[i];
    }
    return 0;
}





