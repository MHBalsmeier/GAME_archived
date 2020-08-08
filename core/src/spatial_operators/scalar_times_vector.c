/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include <stdio.h>

int scalar_times_vector(Scalar_field in_field_0, Vector_field in_field_1, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, lower_index, upper_index;
    double scalar_value, total_volume, upper_volume, lower_volume, upper_weight, lower_weight;
	#pragma omp parallel for private (layer_index, h_index, lower_index, upper_index, scalar_value, total_volume, upper_volume, lower_volume, upper_weight, lower_weight)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            scalar_value = 0.5*(in_field_0[grid -> to_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H] + in_field_0[grid -> from_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]);
        }
        else
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            upper_volume = grid -> volume_ratios[2*upper_index + 1]*grid -> volume[upper_index];
            lower_volume = grid -> volume_ratios[2*lower_index + 0]*grid -> volume[lower_index];
            total_volume = upper_volume + lower_volume;
            upper_weight = upper_volume/total_volume;
            lower_weight = lower_volume/total_volume;
            scalar_value = upper_weight*in_field_0[upper_index] + lower_weight*in_field_0[lower_index];
        }
        out_field[i] = scalar_value*in_field_1[i];
    }
    #pragma omp parallel for private(scalar_value)
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        scalar_value = in_field_0[i] + 0.5*(in_field_0[i] - in_field_0[i + NO_OF_SCALARS_H]);
        out_field[i] = scalar_value*in_field_1[i];
    }
    #pragma omp parallel for private(layer_index, h_index, upper_index, scalar_value)
    for (int i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
        scalar_value = in_field_0[upper_index] + 0.5*(in_field_0[upper_index] - in_field_0[upper_index - NO_OF_SCALARS_H]);
        out_field[i] = scalar_value*in_field_1[i];
    }
    return 0;
}

int scalar_times_vector_for_advection(Scalar_field in_field_0, Vector_field in_field_1, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, lower_index, upper_index;
    double scalar_value, total_volume, upper_volume, lower_volume, upper_weight, lower_weight;
	#pragma omp parallel for private (layer_index, h_index, lower_index, upper_index, scalar_value, total_volume, upper_volume, lower_volume, upper_weight, lower_weight)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            scalar_value = 0.5*(in_field_0[grid -> to_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H] + in_field_0[grid -> from_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]);
        }
        else
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            upper_volume = grid -> volume_ratios[2*upper_index + 1]*grid -> volume[upper_index];
            lower_volume = grid -> volume_ratios[2*lower_index + 0]*grid -> volume[lower_index];
            total_volume = upper_volume + lower_volume;
            upper_weight = upper_volume/total_volume;
            lower_weight = lower_volume/total_volume;
            scalar_value = upper_weight*in_field_0[upper_index] + lower_weight*in_field_0[lower_index];
        }
        out_field[i] = scalar_value*in_field_1[i];
    }
    #pragma omp parallel for private(scalar_value)
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        scalar_value = in_field_0[i] + 0.5*(in_field_0[i] - in_field_0[i + NO_OF_SCALARS_H]);
        out_field[i] = scalar_value*in_field_1[i];
    }
    #pragma omp parallel for private(layer_index, h_index, upper_index, scalar_value)
    for (int i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
        scalar_value = in_field_0[upper_index] + 0.5*(in_field_0[upper_index] - in_field_0[upper_index - NO_OF_SCALARS_H]);
        out_field[i] = scalar_value*in_field_1[i];
    }
    return 0;
}

int scalar_times_vector_scalar_h_v(Scalar_field in_field_h, Scalar_field in_field_v, Vector_field in_field_1, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, lower_index, upper_index;
    double scalar_value, total_volume, upper_volume, lower_volume, upper_weight, lower_weight;
    #pragma omp parallel for private (layer_index, h_index, lower_index, upper_index, scalar_value, total_volume, upper_volume, lower_volume, upper_weight, lower_weight)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            scalar_value = 0.5*(in_field_h[grid -> to_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H] + in_field_h[grid -> from_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]);
        }
        else
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            upper_volume = grid -> volume_ratios[2*upper_index + 1]*grid -> volume[upper_index];
            lower_volume = grid -> volume_ratios[2*lower_index + 0]*grid -> volume[lower_index];
            total_volume = upper_volume + lower_volume;
            upper_weight = upper_volume/total_volume;
            lower_weight = lower_volume/total_volume;
            scalar_value = upper_weight*in_field_v[upper_index] + lower_weight*in_field_v[lower_index];
        }
        out_field[i] = scalar_value*in_field_1[i];
    }
    #pragma omp parallel for private(scalar_value)
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        scalar_value = in_field_v[i] + 0.5*(in_field_v[i] - in_field_v[i + NO_OF_SCALARS_H]);
        out_field[i] = scalar_value*in_field_1[i];
    }
    #pragma omp parallel for private(layer_index, h_index, upper_index, scalar_value)
    for (int i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
        scalar_value = in_field_v[upper_index] + 0.5*(in_field_v[upper_index] - in_field_v[upper_index - NO_OF_SCALARS_H]);
        out_field[i] = scalar_value*in_field_1[i];
    }
    return 0;
}

int scalar_times_vector_vector_h_v(Scalar_field in_field_0, Vector_field in_field_1, Vector_field in_field_2, Vector_field out_field, Grid *grid)
{
	// in_field_1 is horizontal vector field, in_field_2 is vertical vector_field
    int layer_index, h_index, lower_index, upper_index;
    double scalar_value, total_volume, upper_volume, lower_volume, upper_weight, lower_weight;
    #pragma omp parallel for private(layer_index, h_index, lower_index, upper_index, scalar_value, total_volume, upper_volume, lower_volume, upper_weight, lower_weight)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            scalar_value = 0.5*(in_field_0[grid -> to_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H] + in_field_0[grid -> from_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H]);
        }
        else
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            upper_volume = grid -> volume_ratios[2*upper_index + 1]*grid -> volume[upper_index];
            lower_volume = grid -> volume_ratios[2*lower_index + 0]*grid -> volume[lower_index];
            total_volume = upper_volume + lower_volume;
            upper_weight = upper_volume/total_volume;
            lower_weight = lower_volume/total_volume;
            scalar_value = upper_weight*in_field_0[upper_index] + lower_weight*in_field_0[lower_index];
        }
        out_field[i] = scalar_value*in_field_2[i];
    }
    #pragma omp parallel for private(scalar_value)
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        scalar_value = in_field_0[i] + 0.5*(in_field_0[i] - in_field_0[i + NO_OF_SCALARS_H]);
        out_field[i] = scalar_value*in_field_2[i];
    }
    #pragma omp parallel for private(layer_index, h_index, upper_index, scalar_value)
    for (int i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
        scalar_value = in_field_0[upper_index] + 0.5*(in_field_0[upper_index] - in_field_0[upper_index - NO_OF_SCALARS_H]);
        out_field[i] = scalar_value*in_field_2[i];
    }
    return 0;
}






