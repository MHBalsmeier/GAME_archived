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

int calc_kinetic_energy(double latitude_scalar[], double longitude_scalar[], int e_kin_indices[], double e_kin_weights[], double volume[], int adjacent_vector_indices_dual_h[], double alpha_1, int to_index[], int from_index[], double area_dual_pre[], double area[], double z_scalar[], double z_vector[], int adjacent_vector_indices_h[], double latitude_vector[], double longitude_vector[], double latitude_scalar_dual[], double longitude_scalar_dual[], int to_index_dual[], int from_index_dual[], double z_vector_dual[])
{
	double alpha_2 = 1 - alpha_1;
	double z_value_0, z_value_1, z_value_2, z_triangle_mean, partial_volume, kite_area_0, kite_area_1, kite_area, triangle_area_0, triangle_area_1, triangle_area_0_for_kite, triangle_area_1_for_kite, triangle_area_for_kite, triangle_face_0, triangle_face_1, cell_base_area, base_area;
	int second_index, index_found, neighbouring_cell_h_index, vertex_index, layer_index, h_index, number_of_edges, check_result;
	int indices_found_vector[6];
	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
	{
		layer_index = i/NUMBER_OF_SCALARS_H;
		h_index = i - layer_index*NUMBER_OF_SCALARS_H;
		number_of_edges = 6;
		if (h_index < NUMBER_OF_PENTAGONS)
			number_of_edges = 5;
		cell_base_area = area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
		base_area = cell_base_area*pow((RADIUS + z_scalar[i])/(RADIUS + z_vector[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER]), 2);
		for (int j = 0; j < 6; ++j)
		{
			if (j < 5 || h_index >= NUMBER_OF_PENTAGONS)
			{
				e_kin_indices[14*i + j] = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j];
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
				// The Ringler formulation of the kinetic energy is completed here. The rest is Gassmann in order to avoid Hollingsworth instability.
				calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + j]], longitude_vector[adjacent_vector_indices_h[6*h_index + j]], latitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area_0);
				second_index = -1;
				for (int k = 0; k < number_of_edges; ++k)
				{
					if (k != j && (to_index_dual[adjacent_vector_indices_h[6*h_index + k]] == from_index_dual[adjacent_vector_indices_h[6*h_index + j]] || from_index_dual[adjacent_vector_indices_h[6*h_index + k]] == from_index_dual[adjacent_vector_indices_h[6*h_index + j]]))
						second_index = k;
				}
				calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + second_index]], longitude_vector[adjacent_vector_indices_h[6*h_index + second_index]], latitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area_1);
				kite_area_0 = triangle_area_0 + triangle_area_1;
				calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + j]], longitude_vector[adjacent_vector_indices_h[6*h_index + j]], latitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area_0);
				second_index = -1;
				for (int k = 0; k < number_of_edges; ++k)
				{
					if (k != j && (to_index_dual[adjacent_vector_indices_h[6*h_index + k]] == to_index_dual[adjacent_vector_indices_h[6*h_index + j]] || from_index_dual[adjacent_vector_indices_h[6*h_index + k]] == to_index_dual[adjacent_vector_indices_h[6*h_index + j]]))
						second_index = k;
				}
				calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + second_index]], longitude_vector[adjacent_vector_indices_h[6*h_index + second_index]], latitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area_1);
				kite_area_1 = triangle_area_0 + triangle_area_1;
				// mapping kite_area_0 from the unit sphere to the respective model layer height
				kite_area_0 = kite_area_0*pow(RADIUS + z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], 2);
				// mapping kite_area_1 from the unit sphere to the respective model layer height
				kite_area_1 = kite_area_1*pow(RADIUS + z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], 2);
				// finding the neighbouring cell
				neighbouring_cell_h_index = from_index[adjacent_vector_indices_h[6*h_index + j]];
				if (neighbouring_cell_h_index == h_index)
					neighbouring_cell_h_index = to_index[adjacent_vector_indices_h[6*h_index + j]];
				// calculating triangle_area_0 on the unit sphere
				calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_scalar[neighbouring_cell_h_index], longitude_scalar[neighbouring_cell_h_index], latitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area_0_for_kite);
				// calculating triangle_area_1 on the unit sphere
				calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_scalar[neighbouring_cell_h_index], longitude_scalar[neighbouring_cell_h_index], latitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area_1_for_kite);
				// mapping triangle_area_0 from the unit sphere to the model layer height
				triangle_area_0_for_kite = triangle_area_0_for_kite*pow(RADIUS + z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], 2);
				// mapping triangle_area_1 from the unit sphere to the model layer height
				triangle_area_1_for_kite = triangle_area_1_for_kite*pow(RADIUS + z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], 2);
				e_kin_weights[14*i + j] = alpha_1*(triangle_face_0 + triangle_face_1);
				e_kin_weights[14*i + j] += alpha_2*(kite_area_0*triangle_area_0_for_kite/area_dual_pre[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + from_index_dual[adjacent_vector_indices_h[6*h_index + j]]] + kite_area_1*triangle_area_1_for_kite/area_dual_pre[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + to_index_dual[adjacent_vector_indices_h[6*h_index + j]]]);
				e_kin_weights[14*i + j] = e_kin_weights[14*i + j]/base_area;
			}
			else
			{
				e_kin_indices[14*i + j] = 0;
				e_kin_weights[14*i + j] = 0;
			}
		}
		// Now the 6 (5 for pentagonal cells) vector components are taken into account that come into play through the Gassmann formulation of the kinetic energy.
		for (int j = 0; j < 6; ++j)
			indices_found_vector[j] = -2;
		for (int j = 6; j < 12; ++j)
		{
			// the vertex_index is the dual scalar index to which the new vector component belongs
			vertex_index = from_index_dual[adjacent_vector_indices_h[6*h_index + j - 6]];
			// each vertex_index may only appear once per cell
			in_bool_calculator(indices_found_vector, 6, vertex_index, &check_result);
			if (check_result == 1)
			{
				vertex_index = to_index_dual[adjacent_vector_indices_h[6*h_index + j - 6]];
			}
			indices_found_vector[j - 6] = vertex_index;
			// this if is necessary for the pentagonal cells
			if (j < 11 || h_index >= NUMBER_OF_PENTAGONS)
			{
				// the second_index is the second interesection of the triangle with the cell apart from adjacent_vector_indices_h[6*h_index + j - 6]
				second_index = -1;
				for (int k = 0; k < number_of_edges; ++k)
				{
					if (k != j - 6 && (to_index_dual[adjacent_vector_indices_h[6*h_index + k]] == vertex_index || from_index_dual[adjacent_vector_indices_h[6*h_index + k]] == vertex_index))
						second_index = k;
				}
				index_found = 0;
				if (adjacent_vector_indices_dual_h[3*vertex_index + 0] == adjacent_vector_indices_h[6*h_index + j - 6] || adjacent_vector_indices_dual_h[3*vertex_index + 0] == adjacent_vector_indices_h[6*h_index + second_index])
				{
					index_found = 1;
				}
				if (adjacent_vector_indices_dual_h[3*vertex_index + 1] == adjacent_vector_indices_h[6*h_index + j - 6] || adjacent_vector_indices_dual_h[3*vertex_index + 1] == adjacent_vector_indices_h[6*h_index + second_index])
				{
					index_found = 2;
				}
				// the respective vector index is now found
				e_kin_indices[14*i + j] = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + adjacent_vector_indices_dual_h[3*vertex_index + index_found];
				// now the weight must be computed
				// triangle_area_0 is the first part of the kite
				calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + j - 6]], longitude_vector[adjacent_vector_indices_h[6*h_index + j - 6]], latitude_scalar_dual[vertex_index], longitude_scalar_dual[vertex_index], &triangle_area_0);
				// triangle_area_1 is the second part of the kite
				calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + second_index]], longitude_vector[adjacent_vector_indices_h[6*h_index + second_index]], latitude_scalar_dual[vertex_index], longitude_scalar_dual[vertex_index], &triangle_area_1);
				// the kite_area is the sum of triangle_area_0 and triangle_area_1
				kite_area = triangle_area_0 + triangle_area_1;
				// mapping the kite_area from the unit sphere to the respective model layer height
				kite_area = kite_area*pow(RADIUS + z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + vertex_index], 2);
				// calculating triangle_area_for_kite on the uni sphere
				calc_triangle_face(latitude_scalar[to_index[adjacent_vector_indices_dual_h[3*vertex_index + index_found]]], longitude_scalar[to_index[adjacent_vector_indices_dual_h[3*vertex_index + index_found]]], latitude_scalar[from_index[adjacent_vector_indices_dual_h[3*vertex_index + index_found]]], longitude_scalar[from_index[adjacent_vector_indices_dual_h[3*vertex_index + index_found]]], latitude_scalar_dual[vertex_index], longitude_scalar_dual[vertex_index], &triangle_area_for_kite);
				// rescaling triangle_area_for_kite to the respective model layer height
				triangle_area_for_kite = triangle_area_for_kite*pow(RADIUS + z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + vertex_index], 2);
				// calculating the final result for these components
				e_kin_weights[14*i + j] = alpha_2*kite_area/base_area*triangle_area_for_kite/area_dual_pre[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + vertex_index];
			}
			else
			{
				e_kin_indices[14*i + j] = 0;
				e_kin_weights[14*i + j] = 0;    	
			}
		}
		e_kin_indices[14*i + 12] = h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER;
		partial_volume = find_volume(base_area, RADIUS + z_scalar[i], RADIUS + z_vector[e_kin_indices[14*i + 12]]);
		e_kin_weights[14*i + 12] = 0.5*partial_volume/volume[i];
		e_kin_indices[14*i + 13] = h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER;
		partial_volume = find_volume(cell_base_area, RADIUS + z_vector[e_kin_indices[14*i + 13]], RADIUS + z_scalar[i]);
		e_kin_weights[14*i + 13] = 0.5*partial_volume/volume[i];
	}
	double e_kin_weights_sum;
	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
	{
		e_kin_weights_sum = 0;
		for (int j = 0; j < 14; ++j)
		{
			if (e_kin_indices[14*i + j] < 0 || e_kin_indices[14*i + j] >= NUMBER_OF_VECTORS)
			{
				printf("Error in e_kin_indices, position 0.\n");
				exit(1);
			}
			e_kin_weights_sum += e_kin_weights[14*i + j];
			if (e_kin_weights[14*i + j] < 0)
			{
				printf("Error in e_kin_weights, position 0.\n");
				exit(1);
			}
		}    	
		if (fabs(e_kin_weights_sum - 1.5) > 0.015)
		{
			printf("Error in e_kin_weights, position 1. Weights sum is %lf.\n", e_kin_weights_sum);
			for (int j = 0; j < 14; ++j)
			{
				printf("%lf\n", e_kin_weights[14*i + j]);
			}
			printf("%d\n", i);
			exit(1);
		}
	}
	return 0;
}



