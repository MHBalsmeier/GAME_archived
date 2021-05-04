/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
In this file, algebraic multiplications of fields are collected.
*/

#include "../enum_and_typedefs.h"
#include "spatial_operators.h"
#include "../thermodynamics/thermodynamics.h"
#include <stdio.h>

int scalar_times_scalar(Scalar_field in_field_0, Scalar_field in_field_1, Scalar_field out_field)
{
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	out_field[i] = in_field_0[i]*in_field_1[i];
	}
    return 0;
}

int vector_times_vector(Vector_field in_field_0, Vector_field in_field_1, Vector_field out_field)
{
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	out_field[i] = in_field_0[i]*in_field_1[i];
	}
    return 0;
}

int scalar_times_vector(Scalar_field scalar_field, Vector_field vector_field, Vector_field out_field, Grid *grid)
{
	/*
	This function multiplies the vector field vector_field by the scalar field scalar_field.
	*/
    scalar_times_vector_scalar_h(scalar_field, vector_field, out_field, grid);
    scalar_times_vector_scalar_v(scalar_field, vector_field, out_field, grid);
    return 0;
}

int scalar_times_vector_scalar_h(Scalar_field in_field_h, Vector_field vector_field, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, vector_index;
    double scalar_value;
    #pragma omp parallel for private (layer_index, h_index, vector_index, scalar_value)
    for (int i = 0; i < NO_OF_H_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_H;
        h_index = i - layer_index*NO_OF_VECTORS_H;
        vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
        scalar_value
        = 0.5*(
        in_field_h[grid -> to_index[h_index] + layer_index*NO_OF_SCALARS_H]
        + in_field_h[grid -> from_index[h_index] + layer_index*NO_OF_SCALARS_H]);
    	out_field[vector_index] = scalar_value*vector_field[vector_index];
    }
    return 0;
}

int scalar_times_vector_scalar_v(Scalar_field in_field_v, Vector_field vector_field, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, lower_index, upper_index;
    double scalar_value;
    #pragma omp parallel for private (layer_index, h_index, lower_index, upper_index, scalar_value)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index < NO_OF_SCALARS_H)
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            scalar_value = 0.5*(
            in_field_v[upper_index]
            + in_field_v[lower_index]);
        	out_field[i] = scalar_value*vector_field[i];
        }
    }
    // linear extrapolation to the TOA
    #pragma omp parallel for private(scalar_value)
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
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
    for (int i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
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








