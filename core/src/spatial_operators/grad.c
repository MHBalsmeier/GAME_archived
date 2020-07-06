/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "../diagnostics/diagnostics.h"
#include <stdio.h>

int grad(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
    int layer_index, h_index, lower_index, upper_index, retval;
    double vertical_gradient, dzdx;
    for (int i = NO_OF_VECTORS_V; i < NO_OF_VECTORS - NO_OF_VECTORS_V; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_VECTORS_V)
            out_field[i] = (in_field[grid -> to_index[h_index - NO_OF_VECTORS_V] + layer_index*NO_OF_SCALARS_H] - in_field[grid -> from_index[h_index - NO_OF_VECTORS_V] + layer_index*NO_OF_SCALARS_H])/grid -> normal_distance[i];
        else
        {
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            out_field[i] = (in_field[upper_index] - in_field[lower_index])/grid -> normal_distance[i];
        }
    }
    for (int i = NO_OF_VECTORS_V; i < NO_OF_VECTORS - NO_OF_VECTORS_V; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_VECTORS_V && layer_index >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
        {
        	dzdx = (grid -> z_scalar[grid -> to_index[h_index - NO_OF_VECTORS_V] + layer_index*NO_OF_SCALARS_H] - grid -> z_scalar[grid -> from_index[h_index - NO_OF_VECTORS_V] + layer_index*NO_OF_SCALARS_H])/grid -> normal_distance[i];
        	retval = recov_hor_ver_pri(out_field, layer_index, h_index - NO_OF_VECTORS_V, &vertical_gradient, grid);
        	if (retval != 0)
        		printf("Error in recov_hor_ver_pri called from grad.\n");
            out_field[i] = out_field[i] - dzdx*vertical_gradient;
        }
    }
    for (int i = 0; i < NO_OF_VECTORS_V; ++i)
        out_field[i] = out_field[i + NO_OF_VECTORS_PER_LAYER];
    for (int i = NO_OF_VECTORS - NO_OF_VECTORS_V; i < NO_OF_VECTORS; ++i)
        out_field[i] = out_field[i - NO_OF_VECTORS_PER_LAYER];
    return 0;
}

int grad_v_vector_column_to_vector_points(double vertical_velocity[], double grad_vector[], int h_index, Grid *grid)
{
	for (int i = 0; i < NO_OF_LAYERS - 1; ++i)
	{
		if (i == 0)
			grad_vector[i] = -vertical_velocity[1]/(grid -> z_vector[h_index + i*NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[h_index + (i + 2)*NO_OF_VECTORS_PER_LAYER]);
		else if (i == NO_OF_LAYERS - 2)
			grad_vector[i] = vertical_velocity[i - 1]/(grid -> z_vector[h_index + i*NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[h_index + (i + 2)*NO_OF_VECTORS_PER_LAYER]);
		else
			grad_vector[i] = (vertical_velocity[i - 1] - vertical_velocity[i + 1])/(grid -> z_vector[h_index + i*NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[h_index + (i + 2)*NO_OF_VECTORS_PER_LAYER]);
	}
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

int grad_v_scalar_column_to_scalar_points(double scalar_property[], double grad_at_scalars[], int h_index, Grid *grid)
{
	for (int i = 0; i < NO_OF_LAYERS; ++i)
	{
		if (i == 0)
		{
			grad_at_scalars[i] = (scalar_property[i] - scalar_property[i + 1])/(grid -> z_scalar[h_index + i*NO_OF_SCALARS_H] - grid -> z_scalar[h_index + (i + 1)*NO_OF_SCALARS_H]);
		}
		else if (i == NO_OF_LAYERS - 1)
		{
			grad_at_scalars[i] = (scalar_property[i - 1] - scalar_property[i])/(grid -> z_scalar[h_index + (i - 1)*NO_OF_SCALARS_H] - grid -> z_scalar[h_index + i*NO_OF_SCALARS_H]);
		}
		else
		{
			grad_at_scalars[i] = (scalar_property[i - 1] - scalar_property[i + 1])/(grid -> z_scalar[h_index + (i - 1)*NO_OF_SCALARS_H] - grid -> z_scalar[h_index + (i + 1)*NO_OF_SCALARS_H]);
		}
	}
	return 0;
}

int grad_v_vector_column_to_scalar_points(double vertical_velocity[], double grad_vector[], int h_index, Grid *grid)
{
	for (int i = 0; i < NO_OF_LAYERS; ++i)
	{
		if (i == 0)
			grad_vector[i] = -vertical_velocity[0]/(grid -> z_vector[h_index + i*NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[h_index + (i + 1)*NO_OF_VECTORS_PER_LAYER]);
		else if (i == NO_OF_LAYERS - 1)
			grad_vector[i] = vertical_velocity[i - 1]/(grid -> z_vector[h_index + i*NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[h_index + (i + 1)*NO_OF_VECTORS_PER_LAYER]);
		else
			grad_vector[i] = (vertical_velocity[i - 1] - vertical_velocity[i])/(grid -> z_vector[h_index + i*NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[h_index + (i + 1)*NO_OF_VECTORS_PER_LAYER]);
	}
	return 0;
}

















