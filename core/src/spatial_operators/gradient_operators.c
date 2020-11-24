/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "../diagnostics/diagnostics.h"
#include "spatial_operators.h"
#include <stdio.h>

int grad(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, lower_index, upper_index;
    double vertical_gradient;
	#pragma omp parallel for private(layer_index, h_index, lower_index, upper_index, vertical_gradient)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            out_field[i] = (in_field[grid -> to_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H] - in_field[grid -> from_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H])/grid -> normal_distance[i];
        }
        else
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            out_field[i] = (in_field[upper_index] - in_field[lower_index])/grid -> normal_distance[i];
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        out_field[i] = out_field[i + NO_OF_VECTORS_PER_LAYER] + out_field[i + NO_OF_VECTORS_PER_LAYER] - out_field[i + 2*NO_OF_VECTORS_PER_LAYER];
    }
    #pragma omp parallel for
    for (int i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
    {
        out_field[i] = out_field[i - NO_OF_VECTORS_PER_LAYER] + out_field[i - NO_OF_VECTORS_PER_LAYER] - out_field[i - 2*NO_OF_VECTORS_PER_LAYER];
    }
	#pragma omp parallel for private(layer_index, h_index, vertical_gradient)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H && layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
        {
        	recov_hor_ver_pri(out_field, layer_index, h_index - NO_OF_SCALARS_H, &vertical_gradient, grid);
            out_field[i] = out_field[i] - grid -> slope[i]*vertical_gradient;
        }
    }
    return 0;
}

int scalar_times_grad(Scalar_field in_field_for_prefactor, Scalar_field in_field_for_grad, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, lower_index, upper_index;
    double vertical_component;
	#pragma omp parallel for private(layer_index, h_index, lower_index, upper_index, vertical_component)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H)
        {
            out_field[i] = (in_field_for_grad[grid -> to_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H] - in_field_for_grad[grid -> from_index[h_index - NO_OF_SCALARS_H] + layer_index*NO_OF_SCALARS_H])/grid -> normal_distance[i];
        }
        else
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            out_field[i] = (in_field_for_grad[upper_index] - in_field_for_grad[lower_index])/grid -> normal_distance[i];
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        out_field[i] = out_field[i + NO_OF_VECTORS_PER_LAYER] + out_field[i + NO_OF_VECTORS_PER_LAYER] - out_field[i + 2*NO_OF_VECTORS_PER_LAYER];
    }
    #pragma omp parallel for
    for (int i = NO_OF_VECTORS - NO_OF_SCALARS_H; i < NO_OF_VECTORS; ++i)
    {
        out_field[i] = out_field[i - NO_OF_VECTORS_PER_LAYER] + out_field[i - NO_OF_VECTORS_PER_LAYER] - out_field[i - 2*NO_OF_VECTORS_PER_LAYER];
    }
	#pragma omp parallel for private(layer_index, h_index, vertical_component)
    for (int i = NO_OF_SCALARS_H; i < NO_OF_VECTORS - NO_OF_SCALARS_H; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_SCALARS_H && layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
        {
        	recov_hor_ver_pri(out_field, layer_index, h_index - NO_OF_SCALARS_H, &vertical_component, grid);
            out_field[i] = out_field[i] - grid -> slope[i]*vertical_component;
        }
    }
	scalar_times_vector(in_field_for_prefactor, out_field, out_field, grid, 0);
    return 0;
}

int grad_v_scalar_column(double scalar_property[], double grad_vector[], int h_index, Grid *grid)
{
	for (int i = 0; i < NO_OF_LAYERS - 1; ++i)
	{
		grad_vector[i] = (scalar_property[i] - scalar_property[i + 1])/(grid -> z_scalar[h_index + i*NO_OF_SCALARS_H] - grid -> z_scalar[h_index + (i + 1)*NO_OF_SCALARS_H]);
	}
	return 0;
}










