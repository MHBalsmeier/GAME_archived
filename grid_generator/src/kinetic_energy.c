/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
In this file, the kinetic energy indices and weights are computed.
*/

#include "grid_generator.h"
#include "enum.h"
#include <stdlib.h>
#include <stdio.h>
#include "geos95.h"

int calc_kinetic_energy(double latitude_scalar[], double longitude_scalar[], int e_kin_indices[], double e_kin_weights[], double volume[], int adjacent_vector_indices_dual_h[], int to_index[], int from_index[], double area_dual_pre[], double area[], double z_scalar[], double z_vector[], int adjacent_vector_indices_h[], double latitude_vector[], double longitude_vector[], double latitude_scalar_dual[], double longitude_scalar_dual[], int to_index_dual[], int from_index_dual[], double z_vector_dual[])
{
	int layer_index, h_index;
	double z_value_0, z_value_1, z_value_2, z_triangle_mean, triangle_face_0, triangle_face_1, cell_base_area, base_area;
	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
	{
		layer_index = i/NUMBER_OF_SCALARS_H;
		h_index = i - layer_index*NUMBER_OF_SCALARS_H;
		cell_base_area = area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
		base_area = cell_base_area*pow((RADIUS + z_scalar[i])/(RADIUS + z_vector[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER]), 2);
		for (int j = 0; j < 6; ++j)
		{
			if (j < 5 || h_index >= NUMBER_OF_PENTAGONS)
			{
				e_kin_indices[6*i + j] = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j];
				calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + j]], longitude_vector[adjacent_vector_indices_h[6*h_index + j]], latitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_face_0);
				z_value_0 = z_scalar[i];
				z_value_1 = z_vector[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				z_value_2 = z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + from_index_dual[adjacent_vector_indices_h[6*h_index + j]]];
				z_triangle_mean = 1.0/3*(z_value_0 + z_value_1 + z_value_2);
				triangle_face_0 = triangle_face_0*pow(RADIUS + z_triangle_mean, 2);
				calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + j]], longitude_vector[adjacent_vector_indices_h[6*h_index + j]], latitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_face_1);
				z_value_0 = z_scalar[i];
				z_value_1 = z_vector[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				z_value_2 = z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + to_index_dual[adjacent_vector_indices_h[6*h_index + j]]];
				z_triangle_mean = 1.0/3*(z_value_0 + z_value_1 + z_value_2);
				triangle_face_1 = triangle_face_1*pow(RADIUS + z_triangle_mean, 2);
				e_kin_weights[6*i + j] = triangle_face_0 + triangle_face_1;
				e_kin_weights[6*i + j] = e_kin_weights[6*i + j]/base_area;
			}
			else
			{
				e_kin_indices[6*i + j] = 0;
				e_kin_weights[6*i + j] = 0;
			}
		}
	}
	double e_kin_weights_sum;
	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
	{
		e_kin_weights_sum = 0;
		for (int j = 0; j < 6; ++j)
		{
			if (e_kin_indices[6*i + j] < 0 || e_kin_indices[6*i + j] >= NUMBER_OF_VECTORS)
			{
				printf("Error in e_kin_indices, position 0.\n");
				exit(1);
			}
			e_kin_weights_sum += e_kin_weights[6*i + j];
			if (e_kin_weights[6*i + j] < 0)
			{
				printf("Error in e_kin_weights, position 0.\n");
				exit(1);
			}
		}    	
		if (fabs(e_kin_weights_sum - 1.0) > 0.01)
		{
			printf("Error in e_kin_weights, position 1. Weights sum is %lf.\n", e_kin_weights_sum);
			exit(1);
		}
	}
	return 0;
}



