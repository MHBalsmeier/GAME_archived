#include "../enum_and_typedefs.h"
#include <stdio.h>
#include "diagnostics.h"

double inner_elementary(double [], double []);

int vertical_contravariant_normalized(Vector_field in_field, int layer_index, int h_index, Grid *grid, double *result)
{
	if (h_index < 0 || h_index >= NUMBER_OF_VECTORS_V)
		return 1;
	double x_component, y_component;
	int retval = recov_ver_0_pri(in_field, layer_index, h_index, &x_component, grid);
	if (retval != 0)
		printf("Error in recov_ver_pri_0 called at position 0 from horizontal_covariant_normalized.\n");
	retval = recov_ver_1_pri(in_field, layer_index, h_index, &y_component, grid);
	if (retval != 0)
		printf("Error in recov_ver_pri_1 called at position 0 from horizontal_covariant_normalized.\n");
	int vector_index = layer_index*NUMBER_OF_VECTORS_PER_LAYER + h_index;
	if (vector_index < 0 || vector_index >= NUMBER_OF_VECTORS)
		return 2;
	double velocity_vector[3];
	velocity_vector[0] = x_component;
	velocity_vector[1] = y_component;
	velocity_vector[2] = in_field[vector_index];
	double unit_vector[3];
	int layer_index_oro = layer_index - (NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS);
	unit_vector[0] = grid -> vertical_contravar_unit[3*(layer_index_oro*NUMBER_OF_VECTORS_V + h_index) + 0];
	unit_vector[1] = grid -> vertical_contravar_unit[3*(layer_index_oro*NUMBER_OF_VECTORS_V + h_index) + 1];
	unit_vector[2] = grid -> vertical_contravar_unit[3*(layer_index_oro*NUMBER_OF_VECTORS_V + h_index) + 2];
	if (unit_vector[2] < -0.00001 || unit_vector[2] > 1.00001)
		return 3;
	*result = inner_elementary(velocity_vector, unit_vector);
	return 0;
}

int horizontal_covariant_normalized(Vector_field in_field, int layer_index, int h_index, Grid *grid, double *result)
{
	double vertical_component;
	int retval = recov_hor_ver_pri(in_field, layer_index, h_index, &vertical_component, grid);
	if (retval != 0)
		printf("Error in recov_hor_ver_pri called at position 0 from horizontal_covariant_normalized.\n");
	int vector_index = layer_index*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + h_index;
	double delta_z = grid -> z_scalar[layer_index*NUMBER_OF_SCALARS_H + grid -> to_index[h_index]] - grid -> z_scalar[layer_index*NUMBER_OF_SCALARS_H + grid -> from_index[h_index]];
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














