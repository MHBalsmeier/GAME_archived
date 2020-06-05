#include "../enum_and_typedefs.h"
#include <stdio.h>

double inner_elementary(double [], double []);

int vertical_contravariant_normalized(Vector_field in_field, int layer_index, int h_index, Grid *grid, double *result)
{
	double x_component, y_component;
	int retval = recov_ver_pri_0(in_field, layer_index, h_index, &x_component, grid);
	if (retval != 0)
		printf("Error in recov_ver_pri_0 called at position 0 from horizontal_covariant_normalized.\n");
	int retval = recov_ver_pri_1(in_field, layer_index, h_index, &y_component, grid);
	if (retval != 0)
		printf("Error in recov_ver_pri_1 called at position 0 from horizontal_covariant_normalized.\n");
	int vector_index = layer_index*NUMBER_OF_VECTORS_PER_LAYER + h_index;
	double velocity_vector[3];
	velocity_vector[0] = x_component;
	velocity_vector[1] = y_component;
	velocity_vector[2] = in_field[vector_index];
	double unit_vector[3];
	unit_vector[0] = 0;
	unit_vector[1] = 0;
	unit_vector[2] = 1;
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
	double delta_z = z_scalar[grid -> to_index[vector_index]] - z_scalar[grid -> from_index[vector_index]];
	double delta_x = grid -> normal_distance[];
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
