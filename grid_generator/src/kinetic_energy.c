/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
In this file, the kinetic energy indices and weights are computed.
*/

#include <stdlib.h>
#include <stdio.h>
#include "geos95.h"
#include "enum.h"
#include "grid_generator.h"

int calc_kinetic_energy_and_related(double e_kin_weights[], double normal_distance[], double volume[], int to_index[], int from_index[], double area[], double z_scalar[], double z_vector[], int adjacent_vector_indices_h[], double volume_ratios[], double recov_primal2dual_weights[])
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
				e_kin_weights[8*i + j] = area[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				e_kin_weights[8*i + j] = e_kin_weights[8*i + j]*normal_distance[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				e_kin_weights[8*i + j] = e_kin_weights[8*i + j]/(4*volume[i]);
				weights_sum += e_kin_weights[8*i + j];
			}
			else
			{
				e_kin_weights[8*i + j] = 0;
			}
		}
		if (fabs(weights_sum - 1) > 0.2)
		{
			printf("Error in e_kin_weights, position 0. Weights sum is %lf, should be closer to one.\n", weights_sum);
			exit(1);
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
		e_kin_weights[8*i + 6] = area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER]*delta_z/(4*volume[i]);
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
		e_kin_weights[8*i + 7] = area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER]*delta_z/(4*volume[i]);
	}
	int e_kin_h_index_0, e_kin_h_index_1;
	double total_volume, upper_volume, lower_volume, upper_volume_0, lower_volume_0, upper_volume_1, lower_volume_1, check_sum;
	#pragma omp parallel for private(layer_index, e_kin_h_index_0, e_kin_h_index_1, total_volume, upper_volume, lower_volume, upper_volume_0, lower_volume_0, upper_volume_1, lower_volume_1, check_sum)
	for (int i = 0; i < NO_OF_DUAL_H_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		if (layer_index == 0 || layer_index == NO_OF_LAYERS)
		{
			recov_primal2dual_weights[2*i + 0] = 0;
			recov_primal2dual_weights[2*i + 1] = 0;
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
					e_kin_h_index_0 = j;
				if (adjacent_vector_indices_h[6*to_index[h_index] + j] == h_index)
					e_kin_h_index_1 = j;
			}
			if (e_kin_h_index_0 == -1 || e_kin_h_index_1 == -1)
			{
				printf("Index error in calculating recov_primal2dual.\n");
				exit(1);
			}
			upper_volume_0 = e_kin_weights[8*((layer_index - 1)*NO_OF_SCALARS_H + from_index[h_index]) + e_kin_h_index_0]*upper_volume_0;
			upper_volume_1 = e_kin_weights[8*((layer_index - 1)*NO_OF_SCALARS_H + to_index[h_index]) + e_kin_h_index_1]*upper_volume_1;
			lower_volume_0 = e_kin_weights[8*(layer_index*NO_OF_SCALARS_H + from_index[h_index]) + e_kin_h_index_0]*lower_volume_0;
			lower_volume_1 = e_kin_weights[8*(layer_index*NO_OF_SCALARS_H + to_index[h_index]) + e_kin_h_index_1]*lower_volume_1;
			upper_volume = upper_volume_0 + upper_volume_1;
			lower_volume = lower_volume_0 + lower_volume_1;
			total_volume = upper_volume + lower_volume;
			recov_primal2dual_weights[2*i + 0] = upper_volume/total_volume;
			recov_primal2dual_weights[2*i + 1] = lower_volume/total_volume;
			check_sum = recov_primal2dual_weights[2*i + 0] + recov_primal2dual_weights[2*i + 1];
			if (fabs(check_sum - 1) > 1e-10)
			{
				printf("Error in calculating recov_primal2dual. Check sum is %lf, should be 1.\n", check_sum);
				exit(1);
			}
		}
	}
	return 0;
}










