/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
In this file, the inner product weights are computed.
*/

#include <stdlib.h>
#include <stdio.h>
#include "geos95.h"
#include "enum.h"
#include "include.h"

int calc_inner_product_and_related(double inner_product_weights[], double normal_distance[], double volume[], int to_index[], int from_index[], double area[], double z_scalar[], double z_vector[], int adjacent_vector_indices_h[], double volume_ratios[], double remap_horpri2hordual_vector_weights[])
{
	int layer_index, h_index;
	double delta_z, weights_sum, partial_volume;
	#pragma omp parallel for private(layer_index, h_index, delta_z, weights_sum, partial_volume)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		weights_sum = 0;
		for (int j = 0; j < 6; ++j)
		{
			if (j < 5 || h_index >= NO_OF_PENTAGONS)
			{
				inner_product_weights[8*i + j] = area[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				inner_product_weights[8*i + j] = inner_product_weights[8*i + j]*normal_distance[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				inner_product_weights[8*i + j] = inner_product_weights[8*i + j]/(2*volume[i]);
				weights_sum += inner_product_weights[8*i + j];
			}
			else
			{
				inner_product_weights[8*i + j] = 0;
			}
		}
		// upper w, only needed only for diagnostics
		partial_volume = find_volume(area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER]*pow((RADIUS + z_scalar[i])/(RADIUS + z_vector[h_index + layer_index*NO_OF_VECTORS_PER_LAYER]), 2), RADIUS + z_scalar[i], RADIUS + z_vector[h_index + layer_index*NO_OF_VECTORS_PER_LAYER]);
		volume_ratios[2*i + 0] = partial_volume/volume[i];
		if (layer_index == 0)
		{
			delta_z = 2*(z_vector[h_index] - z_scalar[i]);
		}
		else
		{
			delta_z = z_scalar[i - NO_OF_SCALARS_H] - z_scalar[i];
		}
		inner_product_weights[8*i + 6] = area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER]*delta_z/(2*volume[i]);
		// lower w, only needed only for diagnostics
		partial_volume = find_volume(area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER], RADIUS + z_vector[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER], RADIUS + z_scalar[i]);
		volume_ratios[2*i + 1] = partial_volume/volume[i];
		if (layer_index == NO_OF_LAYERS - 1)
		{
			delta_z = 2*(z_scalar[i] - z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + h_index]);
		}
		else
		{
			delta_z = z_scalar[i] - z_scalar[i + NO_OF_SCALARS_H];
		}
		inner_product_weights[8*i + 7] = area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER]*delta_z/(2*volume[i]);
	}
	int e_kin_h_index_0, e_kin_h_index_1;
	double total_volume, upper_volume, lower_volume, upper_volume_0, lower_volume_0, upper_volume_1, lower_volume_1, check_sum;
	#pragma omp parallel for private(layer_index, h_index, e_kin_h_index_0, e_kin_h_index_1, total_volume, upper_volume, lower_volume, upper_volume_0, lower_volume_0, upper_volume_1, lower_volume_1, check_sum)
	for (int i = 0; i < NO_OF_DUAL_H_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		if (layer_index == 0 || layer_index == NO_OF_LAYERS)
		{
			remap_horpri2hordual_vector_weights[2*i + 0] = 0;
			remap_horpri2hordual_vector_weights[2*i + 1] = 0;
		}
		else
		{
			upper_volume_0 = volume[(layer_index - 1)*NO_OF_SCALARS_H + from_index[h_index]];
			upper_volume_1 = volume[(layer_index - 1)*NO_OF_SCALARS_H] + to_index[h_index];
			lower_volume_0 = volume[layer_index*NO_OF_SCALARS_H + from_index[h_index]];
			lower_volume_1 = volume[layer_index*NO_OF_SCALARS_H] + to_index[h_index];
			e_kin_h_index_0 = -1;
			e_kin_h_index_1 = -1;
			for (int j = 0; j < 6; ++j)
			{
				if (adjacent_vector_indices_h[6*from_index[h_index] + j] == h_index)
				{
					e_kin_h_index_0 = j;
				}
				if (adjacent_vector_indices_h[6*to_index[h_index] + j] == h_index)
				{
					e_kin_h_index_1 = j;
				}
			}
			if (e_kin_h_index_0 == -1 || e_kin_h_index_1 == -1)
			{
				printf("Index error in calculating remap_horpri2hordual_vector_weights.\n");
				exit(1);
			}
			upper_volume_0 = inner_product_weights[8*((layer_index - 1)*NO_OF_SCALARS_H + from_index[h_index]) + e_kin_h_index_0]*upper_volume_0;
			upper_volume_1 = inner_product_weights[8*((layer_index - 1)*NO_OF_SCALARS_H + to_index[h_index]) + e_kin_h_index_1]*upper_volume_1;
			lower_volume_0 = inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + from_index[h_index]) + e_kin_h_index_0]*lower_volume_0;
			lower_volume_1 = inner_product_weights[8*(layer_index*NO_OF_SCALARS_H + to_index[h_index]) + e_kin_h_index_1]*lower_volume_1;
			upper_volume = upper_volume_0 + upper_volume_1;
			lower_volume = lower_volume_0 + lower_volume_1;
			total_volume = upper_volume + lower_volume;
			remap_horpri2hordual_vector_weights[2*i + 0] = upper_volume/total_volume;
			remap_horpri2hordual_vector_weights[2*i + 1] = lower_volume/total_volume;
			check_sum = remap_horpri2hordual_vector_weights[2*i + 0] + remap_horpri2hordual_vector_weights[2*i + 1];
			if (fabs(check_sum - 1) > 1e-10)
			{
				printf("Error in calculating remap_horpri2hordual_vector_weights. Check sum is %lf, should be 1.\n", check_sum);
				exit(1);
			}
		}
	}
	return 0;
}










