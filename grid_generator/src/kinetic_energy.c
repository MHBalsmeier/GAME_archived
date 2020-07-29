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

int calc_kinetic_energy_and_volume_ratios(double latitude_scalar[], double longitude_scalar[], int e_kin_indices[], double e_kin_weights[], double volume[], int adjacent_vector_indices_dual_h[], int to_index[], int from_index[], double area_dual_pre[], double area[], double z_scalar[], double z_vector[], int adjacent_vector_indices_h[], double latitude_vector[], double longitude_vector[], double latitude_scalar_dual[], double longitude_scalar_dual[], int to_index_dual[], int from_index_dual[], double z_vector_dual[], double volume_ratios[])
{
	int layer_index, h_index;
	double z_value_0, z_value_1, z_value_2, z_triangle_mean, triangle_face_0, triangle_face_1, delta_z, weights_sum, weights_rescale, partial_volume;
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		weights_sum = 0;
		for (int j = 0; j < 6; ++j)
		{
			if (j < 5 || h_index >= NO_OF_PENTAGONS)
			{
				e_kin_indices[8*i + j] = NO_OF_VECTORS_V + layer_index*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j];
				calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + j]], longitude_vector[adjacent_vector_indices_h[6*h_index + j]], latitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_face_0);
				z_value_0 = z_scalar[i];
				z_value_1 = z_vector[NO_OF_VECTORS_V + layer_index*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				z_value_2 = z_vector_dual[NO_OF_VECTORS_H + layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + from_index_dual[adjacent_vector_indices_h[6*h_index + j]]];
				z_triangle_mean = 1.0/3*(z_value_0 + z_value_1 + z_value_2);
				triangle_face_0 = triangle_face_0*pow(RADIUS + z_triangle_mean, 2);
				calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + j]], longitude_vector[adjacent_vector_indices_h[6*h_index + j]], latitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_face_1);
				z_value_0 = z_scalar[i];
				z_value_1 = z_vector[NO_OF_VECTORS_V + layer_index*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				z_value_2 = z_vector_dual[NO_OF_VECTORS_H + layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + to_index_dual[adjacent_vector_indices_h[6*h_index + j]]];
				z_triangle_mean = 1.0/3*(z_value_0 + z_value_1 + z_value_2);
				triangle_face_1 = triangle_face_1*pow(RADIUS + z_triangle_mean, 2);
				delta_z = z_vector_dual[layer_index*NO_OF_DUAL_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]] - z_vector_dual[(layer_index + 1)*NO_OF_DUAL_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				e_kin_weights[8*i + j] = (triangle_face_0 + triangle_face_1)*delta_z;
				e_kin_weights[8*i + j] = e_kin_weights[8*i + j]/volume[i];
				weights_sum += e_kin_weights[8*i + j];
			}
			else
			{
				e_kin_indices[8*i + j] = 0;
				e_kin_weights[8*i + j] = 0;
			}
		}
		weights_rescale = 1/weights_sum;
		if (fabs(weights_rescale - 1) > 0.12)
		{
			printf("Error in e_kin_weights, position 0. Weights rescale is %lf.\n", weights_rescale);
			exit(1);
		}
		for (int j = 0; j < 6; ++j)
		{
			e_kin_weights[8*i + j] = weights_rescale*e_kin_weights[8*i + j];
		}
		// upper w, only needed only for diagnostics
		e_kin_indices[8*i + 6] = h_index + layer_index*NO_OF_VECTORS_PER_LAYER;
		partial_volume = find_volume(area[e_kin_indices[8*i + 6]]*pow((RADIUS + z_scalar[i])/(RADIUS + z_vector[e_kin_indices[8*i + 6]]), 2), RADIUS + z_scalar[i], RADIUS + z_vector[e_kin_indices[8*i + 6]]);
		e_kin_weights[8*i + 6] = 0.5*partial_volume/volume[i];
		volume_ratios[2*i + 0] = partial_volume/volume[i];
		// lower w, only needed only for diagnostics
		e_kin_indices[8*i + 7] = h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER;
		partial_volume = find_volume(area[e_kin_indices[8*i + 7]], RADIUS + z_vector[e_kin_indices[8*i + 7]], RADIUS + z_scalar[i]);
		e_kin_weights[8*i + 7] = 0.5*partial_volume/volume[i];
		volume_ratios[2*i + 1] = partial_volume/volume[i];
		if (fabs(volume_ratios[2*i + 0] + volume_ratios[2*i + 1] - 1) > 1e-10)
		{
			printf("Problem with volume ratios, volume ratio sum is %lf.\n", volume_ratios[2*i + 0] + volume_ratios[2*i + 1]);
			exit(1);
		}
	}
	double e_kin_weights_check_sum;
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		e_kin_weights_check_sum = 0;
		for (int j = 0; j < 8; ++j)
		{
			if (e_kin_indices[8*i + j] < 0 || e_kin_indices[8*i + j] >= NO_OF_VECTORS)
			{
				printf("Error in e_kin_indices, position 0.\n");
				exit(1);
			}
			e_kin_weights_check_sum += e_kin_weights[8*i + j];
			if (e_kin_weights[8*i + j] < 0)
			{
				printf("Error in e_kin_weights, position 1.\n");
				exit(1);
			}
		}    	
		if (fabs(e_kin_weights_check_sum - 1.5) > 1e-10)
		{
			printf("Error in e_kin_weights, position 1. Weights check sum is %lf.\n", e_kin_weights_check_sum);
			exit(1);
		}
	}
	return 0;
}



