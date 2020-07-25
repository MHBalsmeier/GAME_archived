/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include <stdio.h>
#include "../diagnostics.h"

double inner_elementary(double [], double []);

int vertical_contravariant_normalized(Vector_field in_field, int layer_index, int h_index, Grid *grid, double *result)
{
	if (h_index < 0 || h_index >= NO_OF_SCALARS_H)
		return 1;
	int vector_index = layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
	*result = in_field[vector_index];
	return 0;
}

int vertical_contravariant_normalized_h(Vector_field in_field, int layer_index, int h_index, Grid *grid, double *result)
{
	if (h_index < 0 || h_index >= NO_OF_SCALARS_H)
		return 1;
	// int vector_index = layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
	*result = 0;
	return 0;
}

int horizontal_covariant_normalized(Vector_field in_field, int layer_index, int h_index, Grid *grid, double *result)
{
	double vertical_component;
	recov_hor_ver_pri(in_field, layer_index, h_index, &vertical_component, grid);
	int vector_index = layer_index*NO_OF_VECTORS_PER_LAYER + NO_OF_SCALARS_H + h_index;
	double delta_z = grid -> z_scalar[layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]] - grid -> z_scalar[layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]];
	double delta_x = grid -> normal_distance[vector_index];
	double angle = atan(delta_z/delta_x);
	double velocity_vector[3];
	velocity_vector[0] = in_field[vector_index];
	velocity_vector[1] = 0;
	velocity_vector[2] = vertical_component;
	double unit_vector[3];
	unit_vector[0] = cos(angle);
	unit_vector[1] = 0;
	unit_vector[2] = sin(angle);
	*result = inner_elementary(velocity_vector, unit_vector);
	return 0;
}

double inner_elementary(double vector_0[], double vector_1[])
{
	double result = 0;
	for (int i = 0; i < 3; ++i)
		result += vector_0[i]*vector_1[i];
	return result;
}














