/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
In this file, vector reconstruction indices and weights are computed.
*/

#include "grid_generator.h"
#include "enum.h"
#include <stdlib.h>
#include <stdio.h>
#include "geos95.h"

int calc_coriolis_weights(int from_index_dual[], int to_index_dual[], int trsk_modified_curl_indices[], double normal_distance[], double normal_distance_dual[], int to_index[], double area[], double z_scalar[], double latitude_scalar[], double longitude_scalar[], double latitude_vector[], double longitude_vector[], double latitude_scalar_dual[], double longitude_scalar_dual[], double trsk_modified_weights[], int trsk_modified_velocity_indices[], int from_index[], int adjacent_vector_indices_h[], double direction[], double direction_dual[], double ORTH_CRITERION_DEG, double z_vector[], double z_vector_dual[])
{
	/*
	This function implements the modified TRSK scheme proposed by Gassmann (2018). Indices and weights are computed here for the highest layer but remain unchanged elsewhere.
	*/
	int *face_of_cell_indices = malloc(2*sizeof(int));
	int *from_or_to_index = malloc(NO_OF_VECTORS_H*sizeof(int));
	int offset, sign_0, sign_1, no_of_edges, index_offset;
	double check_sum;
	double rescale_for_z_offset_1d = (RADIUS + z_scalar[0])/(RADIUS + z_vector[0]);
	double rescale_for_z_offset_2d = pow(rescale_for_z_offset_1d, 2);
	double sum_of_weights = 0;
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
		/*
		translation from TRSK paper (Thuburn et al., 2009):
		sign_0: t_{e, v_2}
		sign_1: n_{e', i}
		trsk_modified_weights: w
		*/
		double check_sum, triangle_0, triangle_1;
		int vertex_index_candidate_0, vertex_index_candidate_1, counter, check_result, first_index, last_index;
		first_index = -1;
		last_index = -1;
		for (int k = 0; k < 10; ++k)
		{
			if (k == 0 || k == 5)
				offset = 0;
			if (k < 5)
			{
				index_offset = 0;
				sign_0 = -1;
				for (int l = 0; l < NO_OF_VECTORS_H; ++l)
				{
					from_or_to_index[l] = from_index[l];
				}
			}
			else
			{
				index_offset = 5;
				sign_0 = 1;
				for (int l = 0; l < NO_OF_VECTORS_H; ++l)
				{
					from_or_to_index[l] = to_index[l];
				}
			}
			if (adjacent_vector_indices_h[6*from_or_to_index[i] + k - index_offset] == i)
				offset += 1;
			if (offset > 1)
			{
				printf("Problem 1 in TRSK implementation detected.\n");
				exit(1);
			}
			trsk_modified_velocity_indices[10*i + k] = adjacent_vector_indices_h[6*from_or_to_index[i] + k - index_offset + offset];
			if (trsk_modified_velocity_indices[10*i + k] == -1)
			{
				trsk_modified_weights[10*i + k] = 0;
			}
			else
			{
				sign_1 = -1;
				if (from_index[trsk_modified_velocity_indices[10*i + k]] == from_or_to_index[i])
					sign_1 = 1;
				if (from_or_to_index[i] < NO_OF_PENTAGONS)
				{
					no_of_edges = 5;
				}
				else
				{
					no_of_edges = 6;
				}
				int vertex_indices[no_of_edges];
				int edge_indices[no_of_edges];
				int indices_resorted[no_of_edges];
				int vertex_indices_resorted[no_of_edges];
				double latitude_vertices[no_of_edges];
				double longitude_vertices[no_of_edges];
				double latitude_edges[no_of_edges];
				double longitude_edges[no_of_edges];
				double vector_of_areas[no_of_edges];
				for (int l = 0; l < no_of_edges; ++l)
					vertex_indices[l] = -1;
				counter = 0;
				for (int l = 0; l < no_of_edges; ++l)
				{
					vertex_index_candidate_0 = from_index_dual[adjacent_vector_indices_h[6*from_or_to_index[i] + l]];
					vertex_index_candidate_1 = to_index_dual[adjacent_vector_indices_h[6*from_or_to_index[i] + l]];
					check_result = in_bool_calculator(vertex_index_candidate_0, vertex_indices, no_of_edges);						
					if (check_result == 0)
					{
						vertex_indices[counter] = vertex_index_candidate_0;
						latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
						longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
						++counter;
					}
					check_result = in_bool_calculator(vertex_index_candidate_1, vertex_indices, no_of_edges);						
					if (check_result == 0)
					{
						vertex_indices[counter] = vertex_index_candidate_1;
						latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
						longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
						++counter;
					}
				}
				if (counter != no_of_edges)
				{
					printf("Problem 13 in TRSK implementation detected.\n");
					exit(1);
				}
				sort_edge_indices(latitude_vertices, longitude_vertices, no_of_edges, indices_resorted);
				for (int l = 0; l < no_of_edges; ++l)
					vertex_indices_resorted[l] = vertex_indices[indices_resorted[l]];
				for (int l = 0; l < no_of_edges; ++l)
				{
					for (int m = 0; m < no_of_edges; ++m)
					{
						if ((from_index_dual[adjacent_vector_indices_h[6*from_or_to_index[i] + m]] == vertex_indices_resorted[l] && to_index_dual[adjacent_vector_indices_h[6*from_or_to_index[i] + m]] == vertex_indices_resorted[(l + 1)%no_of_edges]) || (to_index_dual[adjacent_vector_indices_h[6*from_or_to_index[i] + m]] == vertex_indices_resorted[l] && from_index_dual[adjacent_vector_indices_h[6*from_or_to_index[i] + m]] == vertex_indices_resorted[(l + 1)%no_of_edges]))
							edge_indices[l] = adjacent_vector_indices_h[6*from_or_to_index[i] + m];
					}
				}
				for (int l = 0; l < no_of_edges; ++l)
				{
					latitude_edges[l] = latitude_vector[edge_indices[l]];
					longitude_edges[l] = longitude_vector[edge_indices[l]];
				}
				check_sum = 0;
				for (int l = 0; l < no_of_edges; ++l)	
				{
					if (l == 0)
						calc_triangle_face(latitude_scalar[from_or_to_index[i]], longitude_scalar[from_or_to_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[no_of_edges - 1], longitude_edges[no_of_edges - 1], &triangle_0);
					else
						calc_triangle_face(latitude_scalar[from_or_to_index[i]], longitude_scalar[from_or_to_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[l - 1], longitude_edges[l - 1], &triangle_0);
					calc_triangle_face(latitude_scalar[from_or_to_index[i]], longitude_scalar[from_or_to_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[l], longitude_edges[l], &triangle_1);
					vector_of_areas[l] = pow(RADIUS + z_vector[NO_OF_SCALARS_H + i], 2)*(triangle_0 + triangle_1);
					check_sum += vector_of_areas[l];
				}
				if (fabs(check_sum/(rescale_for_z_offset_2d*area[from_or_to_index[i]]) - 1) > 1e-10)
				{
					printf("Problem 30 in TRSK implementation detected. %lf\n", check_sum/(rescale_for_z_offset_2d*area[from_or_to_index[i]]));
					exit(1);
				}
				for (int l = 0; l < no_of_edges; ++l)
				{
					if (edge_indices[l] == i)
						last_index = l;
					if (edge_indices[l] == trsk_modified_velocity_indices[10*i + k])
						first_index = (l + 1)%no_of_edges;
				}
				if (k == index_offset + no_of_edges - 1)
					sum_of_weights = 0;
				else
					double_sum_gen(vector_of_areas, no_of_edges, first_index, last_index, &sum_of_weights);
				if (sum_of_weights < 0 || sum_of_weights/(rescale_for_z_offset_2d*area[from_or_to_index[i]]) > 1)
				{
					printf("Problem 34 in TRSK implementation detected.\n");
					exit(1);
				}
				sum_of_weights = sum_of_weights/(rescale_for_z_offset_2d*area[from_or_to_index[i]]);
				trsk_modified_weights[10*i + k] = sign_0*(sum_of_weights - 0.5)*sign_1;
				trsk_modified_weights[10*i + k] = -rescale_for_z_offset_1d*normal_distance_dual[trsk_modified_velocity_indices[10*i + k]]/normal_distance[NO_OF_SCALARS_H + i]*trsk_modified_weights[10*i + k];
			}
		}
		// modification following Gassmann (2018)
		// First off all, the indices need to be resorted.
		// As usual, the from cell is treated first.
		// First of all it needs to be determined wether the cell at hand is pentagonal or hexagonal.
		no_of_edges = 6;
		if (from_index[i] < NO_OF_PENTAGONS)
			no_of_edges = 5;
		int next_vertex_index, next_vertex_index_candidate;
		double direction_change;
		next_vertex_index = to_index_dual[i];
        find_angle_change(direction[i], direction_dual[i], &direction_change);
        if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
            next_vertex_index = from_index_dual[i];
        int indices_used[no_of_edges - 1];
        int indices_used_counter = 0;
		for (int j = 0; j < no_of_edges - 1; ++j)
		{
			indices_used[j] = -1;
		}
		int trsk_modified_velocity_indices_pre[10];
		double trsk_modified_weights_pre[10];
		for (int j = 0; j < 10; ++j)
		{
			trsk_modified_velocity_indices_pre[j] = trsk_modified_velocity_indices[10*i + j];
			trsk_modified_weights_pre[j] = trsk_modified_weights[10*i + j];
		}
		int value_written;
		for (int j = 0; j < no_of_edges - 1; ++j)
		{
			value_written = 0;
            for (int k = 0; k < no_of_edges - 1; ++k)
            {
            	if ((from_index_dual[trsk_modified_velocity_indices_pre[k]] == next_vertex_index || to_index_dual[trsk_modified_velocity_indices_pre[k]] == next_vertex_index) && 0 == in_bool_calculator(k, indices_used, no_of_edges - 1) && value_written == 0)
            	{
					trsk_modified_velocity_indices[10*i + j] = trsk_modified_velocity_indices_pre[k];
					trsk_modified_weights[10*i + j] = trsk_modified_weights_pre[k];
					indices_used[indices_used_counter] = k;
					indices_used_counter++;
					value_written = 1;
				}
			}
			next_vertex_index_candidate = to_index_dual[trsk_modified_velocity_indices[10*i + j]];
            if (next_vertex_index_candidate == next_vertex_index)
            {
                next_vertex_index = from_index_dual[trsk_modified_velocity_indices[10*i + j]];
            }
            else
            	next_vertex_index = next_vertex_index_candidate;
		}
		if (indices_used_counter != no_of_edges - 1)
		{
			printf("Problem 42 in TRSK implementation detected.\n");
			exit(1);
		}
		// Then comes the to cell.
		// First of all it needs to be determined wether the cell at hand is pentagonal or hexagonal.
		no_of_edges = 6;
		if (to_index[i] < NO_OF_PENTAGONS)
			no_of_edges = 5;
		next_vertex_index = from_index_dual[i];
        find_angle_change(direction[i], direction_dual[i], &direction_change);
        if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
            next_vertex_index = to_index_dual[i];
        indices_used_counter = 0;
		for (int j = 0; j < no_of_edges - 1; ++j)
		{
			indices_used[j] = -1;
		}
		for (int j = 0; j < no_of_edges - 1; ++j)
		{
			value_written = 0;
            for (int k = 0; k < no_of_edges - 1; ++k)
            {
            	if ((from_index_dual[trsk_modified_velocity_indices_pre[5 + k]] == next_vertex_index || to_index_dual[trsk_modified_velocity_indices_pre[5 + k]] == next_vertex_index) && 0 == in_bool_calculator(k, indices_used, no_of_edges - 1) && value_written == 0)
            	{
					trsk_modified_velocity_indices[10*i + 5 + j] = trsk_modified_velocity_indices_pre[5 + k];
					trsk_modified_weights[10*i + 5 + j] = trsk_modified_weights_pre[5 + k];
					indices_used[indices_used_counter] = k;
					indices_used_counter++;
					value_written = 1;
				}
			}
			next_vertex_index_candidate = to_index_dual[trsk_modified_velocity_indices[10*i + 5 + j]];
            if (next_vertex_index_candidate == next_vertex_index)
                next_vertex_index = from_index_dual[trsk_modified_velocity_indices[10*i + 5 + j]];
            else
            	next_vertex_index = next_vertex_index_candidate;
		}
		if (indices_used_counter != no_of_edges - 1)
		{
			printf("Problem 43 in TRSK implementation detected.\n");
			exit(1);
		}
		// Now the resorting itself can be executed.
		if (to_index[i] < NO_OF_PENTAGONS)
			trsk_modified_curl_indices[10*i + 0] = trsk_modified_velocity_indices[10*i + 8];
		else
			trsk_modified_curl_indices[10*i + 0] = trsk_modified_velocity_indices[10*i + 9];
		trsk_modified_curl_indices[10*i + 1] = trsk_modified_velocity_indices[10*i + 0];
		if (from_index[i] < NO_OF_PENTAGONS)
		{
			trsk_modified_curl_indices[10*i + 2] = trsk_modified_velocity_indices[10*i + 3];
			trsk_modified_curl_indices[10*i + 3] = trsk_modified_velocity_indices[10*i + 5];
			trsk_modified_curl_indices[10*i + 4] = 0;
			if (trsk_modified_weights[10*i + 4] != 0)
			{
				printf("Problem 40 in TRSK implementation detected.\n");
				exit(1);
			}
		}
		else
		{
			trsk_modified_curl_indices[10*i + 2] = trsk_modified_velocity_indices[10*i + 2];
			trsk_modified_curl_indices[10*i + 3] = trsk_modified_velocity_indices[10*i + 4];
			trsk_modified_curl_indices[10*i + 4] = trsk_modified_velocity_indices[10*i + 5];
		}
		if (from_index[i] < NO_OF_PENTAGONS)
			trsk_modified_curl_indices[10*i + 5] = trsk_modified_velocity_indices[10*i + 3];
		else
			trsk_modified_curl_indices[10*i + 5] = trsk_modified_velocity_indices[10*i + 4];
		trsk_modified_curl_indices[10*i + 6] = trsk_modified_velocity_indices[10*i + 5];
		if (to_index[i] < NO_OF_PENTAGONS)
		{
			trsk_modified_curl_indices[10*i + 7] = trsk_modified_velocity_indices[10*i + 8];
			trsk_modified_curl_indices[10*i + 8] = trsk_modified_velocity_indices[10*i + 0];
			trsk_modified_curl_indices[10*i + 9] = 0;
			if (trsk_modified_weights[10*i + 9] != 0)
			{
				printf("Problem 41 in TRSK implementation detected.\n");
				exit(1);
			}
		}
		else
		{
			trsk_modified_curl_indices[10*i + 7] = trsk_modified_velocity_indices[10*i + 7];
			trsk_modified_curl_indices[10*i + 8] = trsk_modified_velocity_indices[10*i + 9];
			trsk_modified_curl_indices[10*i + 9] = trsk_modified_velocity_indices[10*i + 0];
		}
		for (int j = 0; j < 10; ++j)
		{
			for (int k = j + 1; k < 10; ++k)
			{
				if (trsk_modified_velocity_indices[10*i + j] == trsk_modified_velocity_indices[10*i + k] && (trsk_modified_weights[10*i + j] != 0 && trsk_modified_weights[10*i + k] != 0))
				{
					printf("Problem 29 in TRSK implementation detected.\n");
					exit(1);
				}
			}
		}
    }
    int first_index, second_index, k;
    // This checks Eq. (39) of the first TRSK paper (Thuburn et al., 2009).
	double value_0, value_1;
	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			first_index = trsk_modified_velocity_indices[10*i + j];
			if (first_index != -1)
			{
				value_0 = normal_distance[NO_OF_SCALARS_H + i]/(rescale_for_z_offset_1d*normal_distance_dual[first_index])*trsk_modified_weights[10*i + j];
				second_index = -1;
				for (k = 0; k < 10; ++k)
				{
					if (trsk_modified_velocity_indices[10*first_index + k] == i)
						second_index = 10*first_index + k;
				}
				if (second_index == -1)
				{
					printf("Problem 38 in TRSK implementation detected.\n");
					exit(1);
				}
				value_1 = normal_distance[NO_OF_SCALARS_H + first_index]/(rescale_for_z_offset_1d*normal_distance_dual[i])*trsk_modified_weights[second_index];
				check_sum = value_0 + value_1;
				if (fabs(check_sum) > 1e-10)
				{
					printf("Problem 39 in TRSK implementation detected.%lf\n", check_sum);
					exit(1);
				}
			}
		}
	}
	for (int i = 0; i < 10*NO_OF_VECTORS_H; ++i)
	{
		if (trsk_modified_velocity_indices[i] == -1)
		{
			trsk_modified_velocity_indices[i] = 0;
		}
		if (trsk_modified_curl_indices[i] == -1)
		{
			trsk_modified_curl_indices[i] = 0;
		}
	}
	free(from_or_to_index);
    free(face_of_cell_indices);
	return 0;
}

int set_recov_ver(int adjacent_vector_indices_h[], double direction[], double direction_dual[], double latitude_scalar[], double longitude_scalar[], double latitude_scalar_dual[], double longitude_scalar_dual[], int from_index_dual[], int to_index_dual[], double pent_hex_face_unity_sphere[], double recov_ver_weight[], double ORTH_CRITERION_DEG, double normal_distance[], double z_vector[], double z_vector_dual[], double TOA)
{
	int no_of_edges, sign, h_index, layer_index;
	double triangle_area, check_sum_pre, check_sum, direction_change, delta_z_at_cell, delta_z_at_edge;
	for (int i = 0; i < NO_OF_LEVELS*NO_OF_SCALARS_H; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		delta_z_at_cell = normal_distance[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
		no_of_edges = 6;
		if (h_index < NO_OF_PENTAGONS)
		{
			no_of_edges = 5;
	    }
	    check_sum_pre = 0;
		for (int j = 0; j < no_of_edges; ++j)
		{
			calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], latitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area);
			sign = 1;
            find_angle_change(direction[adjacent_vector_indices_h[6*h_index + j]], direction_dual[adjacent_vector_indices_h[6*h_index + j]], &direction_change);
            if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
                sign = -1;
            if (layer_index == 0)
            {
            	delta_z_at_edge = TOA - z_vector[NO_OF_SCALARS_H + adjacent_vector_indices_h[6*h_index + j]];
            }
            else if (layer_index == NO_OF_LAYERS)
            {
            	delta_z_at_edge = z_vector[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]] - z_vector_dual[NO_OF_DUAL_VECTORS - NO_OF_VECTORS_H + adjacent_vector_indices_h[6*h_index + j]];
            }
            else
            {
            	delta_z_at_edge = z_vector[NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]] - z_vector[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
            }
		    recov_ver_weight[6*i + j] = 2*sign*triangle_area*delta_z_at_edge/(pent_hex_face_unity_sphere[h_index]*delta_z_at_cell);
		    check_sum_pre += fabs(recov_ver_weight[6*i + j]);
		}
		if (fabs(check_sum_pre - 2) > 0.3)
		{
			printf("Problem with recov_ver_weight, check_sum_pre is %lf.\n", check_sum_pre);
			exit(1);
		}
		check_sum = 0;
		for (int j = 0; j < no_of_edges; ++j)
		{
		    recov_ver_weight[6*i + j] = 2/check_sum_pre*recov_ver_weight[6*i + j];
		    check_sum += fabs(recov_ver_weight[6*i + j]);
		}
		if (fabs(check_sum - 2) > 1e-10)
		{
			printf("Problem with recov_ver_weight, check_sum is %lf.\n", check_sum);
			exit(1);
		}
	}
	return 0;
}






