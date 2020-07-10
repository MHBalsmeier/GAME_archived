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

int calc_coriolis_weights(int recov_hor_ver_curl_index[], int from_index_dual[], int to_index_dual[], double recov_hor_ver_curl_weight[], int recov_hor_ver_pri_index[], int trsk_modified_curl_indices[], double normal_distance[], double normal_distance_dual[], int to_index[], double area[], double z_scalar[], double latitude_scalar[], double longitude_scalar[], double latitude_vector[], double longitude_vector[], double latitude_scalar_dual[], double longitude_scalar_dual[], double trsk_modified_weights[], int trsk_modified_velocity_indices[], int from_index[], int adjacent_vector_indices_h[], double direction[], double recov_hor_par_curl_weight[], double direction_dual[], double rel_on_line_dual[], int recov_hor_par_curl_index[], double ORTH_CRITERION_DEG)
{
	/*
	This function implements the modified TRSK scheme proposed by Gassmann (2018).
	*/
	int *face_of_cell_indices = malloc(2*sizeof(int));
	int offset, sign_0, sign_1, retval, sign, number_of_edges;
	double check_sum, direction_change;
	double sum_of_weights = 0;
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            if (j == 0)
            {
                recov_hor_ver_curl_index[2*i + j] = to_index_dual[i];
                recov_hor_ver_curl_weight[2*i + j] = rel_on_line_dual[i];
            }
            else
            {
                recov_hor_ver_curl_index[2*i + j] = from_index_dual[i];
                recov_hor_ver_curl_weight[2*i + j] = 1 - recov_hor_ver_curl_weight[2*i + 0];
            }
            sign = 1;
            recov_hor_par_curl_index[2*i + j] = i + j*NUMBER_OF_DUAL_VECTORS_PER_LAYER;
            find_angle_change(direction[i], direction_dual[i], &direction_change);
            if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
                sign = -1;
            recov_hor_par_curl_weight[2*i + j] = sign*0.5;
        }
		/*
		translation from TRSK paper:
		sign_0: t_{e, v_2}
		sign_1: n_{e', i}
		trsk_modified_weights: w
		*/
		offset = 0;
		double check_sum, triangle_0, triangle_1;
		int vertex_index_candidate_0, vertex_index_candidate_1, counter, check_result, first_index, last_index;
		first_index = -1;
		last_index = -1;
		for (int k = 0; k < 10; ++k)
		{
			if (k < 5)
			{
				sign_0 = -1;
				if (adjacent_vector_indices_h[6*from_index[i] + k] == i)
					offset += 1;
				if (offset > 1)
				{
					printf("Problem 1 in TRSK implementation detected.\n");
					exit(1);
				}
				trsk_modified_velocity_indices[10*i + k] = adjacent_vector_indices_h[6*from_index[i] + k + offset];
				if (trsk_modified_velocity_indices[10*i + k] == -1)
				{
					trsk_modified_velocity_indices[10*i + k] = 0;
					trsk_modified_weights[10*i + k] = 0;
				}
				else
				{
					if (trsk_modified_velocity_indices[10*i + k] < 0 || trsk_modified_velocity_indices[10*i + k] >= NUMBER_OF_VECTORS_H)
					{
						printf("Problem 11 in TRSK implementation detected.\n");
						exit(1);
					}
					sign_1 = -1;
					if (from_index[trsk_modified_velocity_indices[10*i + k]] == from_index[i])
						sign_1 = 1;
					if (from_index[i] < NUMBER_OF_PENTAGONS)
					{
						number_of_edges = 5;
					}
					else
					{
						number_of_edges = 6;
					}
					int vertex_indices[number_of_edges];
					int edge_indices[number_of_edges];
					int indices_resorted[number_of_edges];
					int vertex_indices_resorted[number_of_edges];
					double latitude_vertices[number_of_edges];
					double longitude_vertices[number_of_edges];
					double latitude_edges[number_of_edges];
					double longitude_edges[number_of_edges];
					double vector_of_areas[number_of_edges];
					for (int l = 0; l < number_of_edges; ++l)
						vertex_indices[l] = -1;
					counter = 0;
					for (int l = 0; l < number_of_edges; ++l)
					{
						vertex_index_candidate_0 = from_index_dual[adjacent_vector_indices_h[6*from_index[i] + l]];
						vertex_index_candidate_1 = to_index_dual[adjacent_vector_indices_h[6*from_index[i] + l]];
						retval = in_bool_calculator(vertex_indices, number_of_edges, vertex_index_candidate_0, &check_result);						
						if (check_result == 0)
						{
							vertex_indices[counter] = vertex_index_candidate_0;
							latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
							longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
							++counter;
						}
						retval = in_bool_calculator(vertex_indices, number_of_edges, vertex_index_candidate_1, &check_result);						
						if (check_result == 0)
						{
							vertex_indices[counter] = vertex_index_candidate_1;
							latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
							longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
							++counter;
						}
					}
					if (counter != number_of_edges)
					{
						printf("Problem 13 in TRSK implementation detected.\n");
						exit(1);
					}
					for (int l = 0; l < number_of_edges; ++l)
					{
						if (vertex_indices[l] < 0 || vertex_indices[l] >= NUMBER_OF_DUAL_SCALARS_H)
						{
							printf("Problem 3 in TRSK implementation detected.\n");
							exit(1);
						}
						for (int m = l + 1; m < number_of_edges; ++m)
						{
							if (vertex_indices[l] == vertex_indices[m])
							{
								printf("Problem 17 in TRSK implementation detected.\n");
								exit(1);
							}
						}
					}
					retval = sort_edge_indices(latitude_vertices, longitude_vertices, number_of_edges, indices_resorted);
					for (int l = 0; l < number_of_edges; ++l)
					{
						for (int m = l + 1; m < number_of_edges; ++m)
						{
							if (indices_resorted[l] == indices_resorted[m] || indices_resorted[l] < 0 || indices_resorted[l] > number_of_edges - 1)
							{
								printf("Problem 25 in TRSK implementation detected.\n");
								exit(1);
							}
						}
					}
					for (int l = 0; l < number_of_edges; ++l)
						vertex_indices_resorted[l] = vertex_indices[indices_resorted[l]];
					for (int l = 0; l < number_of_edges; ++l)
					{
						for (int m = 0; m < number_of_edges; ++m)
						{
							if ((from_index_dual[adjacent_vector_indices_h[6*from_index[i] + m]] == vertex_indices_resorted[l] && to_index_dual[adjacent_vector_indices_h[6*from_index[i] + m]] == vertex_indices_resorted[(l + 1)%number_of_edges]) || (to_index_dual[adjacent_vector_indices_h[6*from_index[i] + m]] == vertex_indices_resorted[l] && from_index_dual[adjacent_vector_indices_h[6*from_index[i] + m]] == vertex_indices_resorted[(l + 1)%number_of_edges]))
								edge_indices[l] = adjacent_vector_indices_h[6*from_index[i] + m];
						}
					}
					for (int l = 0; l < number_of_edges; ++l)
					{
						if (edge_indices[l] < 0 || edge_indices[l] >= NUMBER_OF_VECTORS_H)
						{
							printf("Problem 7 in TRSK implementation detected.\n");
							exit(1);
						}
						for (int m = l + 1; m < number_of_edges; ++m)
						{
							if (edge_indices[l] == edge_indices[m])
							{
								printf("Problem 18 in TRSK implementation detected.\n");
								exit(1);
							}
						}
					}
					for (int l = 0; l < number_of_edges; ++l)
					{
						latitude_edges[l] = latitude_vector[edge_indices[l]];
						longitude_edges[l] = longitude_vector[edge_indices[l]];
					}
					check_sum = 0;
					for (int l = 0; l < number_of_edges; ++l)	
					{
						if (l == 0)
							retval = calc_triangle_face(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[number_of_edges - 1], longitude_edges[number_of_edges - 1], &triangle_0);
						else
							retval = calc_triangle_face(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[l - 1], longitude_edges[l - 1], &triangle_0);
						retval = calc_triangle_face(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[l], longitude_edges[l], &triangle_1);
						vector_of_areas[l] = pow(RADIUS + z_scalar[from_index[i]], 2)*(triangle_0 + triangle_1);
						check_sum += vector_of_areas[l];
					}
					if (fabs(check_sum/area[from_index[i]] - 1) > 0.001)
					{
						printf("Problem 30 in TRSK implementation detected.\n");
						exit(1);
					}
					for (int l = 0; l < number_of_edges; ++l)
					{
						if (edge_indices[l] == i)
							last_index = l;
						if (edge_indices[l] == trsk_modified_velocity_indices[10*i + k])
							first_index = (l + 1)%number_of_edges;
					}
					if (k == number_of_edges - 1)
						sum_of_weights = 0;
					else
						double_sum_gen(vector_of_areas, number_of_edges, first_index, last_index, &sum_of_weights);
					if (sum_of_weights < 0 || sum_of_weights/area[from_index[i]] > 1)
					{
						printf("Problem 34 in TRSK implementation detected\n");
						exit(1);
					}
					sum_of_weights = sum_of_weights/area[from_index[i]];
					trsk_modified_weights[10*i + k] = sign_0*(sum_of_weights - 0.5)*sign_1;
				}
			}
			else
			{
				if (k == 5)
					offset = 0;
				sign_0 = 1;
				if (adjacent_vector_indices_h[6*to_index[i] + k - 5] == i)
					offset += 1;
				if (offset > 1)
				{
					printf("Problem 2 in TRSK implementation detected.\n");
					exit(1);
				}
				trsk_modified_velocity_indices[10*i + k] = adjacent_vector_indices_h[6*to_index[i] + k - 5 + offset];
				if (trsk_modified_velocity_indices[10*i + k] == -1)
				{
					trsk_modified_velocity_indices[10*i + k] = 0;
					trsk_modified_weights[10*i + k] = 0;
				}
				else
				{
					if (trsk_modified_velocity_indices[10*i + k] < 0 || trsk_modified_velocity_indices[10*i + k] >= NUMBER_OF_VECTORS_H)
					{
						printf("Problem 12 in TRSK implementation detected.\n");
						exit(1);
					}
					sign_1 = -1;
					if (from_index[trsk_modified_velocity_indices[10*i + k]] == to_index[i])
						sign_1 = 1;
					if (to_index[i] < NUMBER_OF_PENTAGONS)
					{
						number_of_edges = 5;
					}
					else
					{
						number_of_edges = 6;
					}
					int vertex_indices[number_of_edges];
					int edge_indices[number_of_edges];
					int indices_resorted[number_of_edges];
					int vertex_indices_resorted[number_of_edges];
					double latitude_vertices[number_of_edges];
					double longitude_vertices[number_of_edges];
					double latitude_edges[number_of_edges];
					double longitude_edges[number_of_edges];
					double vector_of_areas[number_of_edges];
					for (int l = 0; l < number_of_edges; ++l)
						vertex_indices[l] = -1;
					counter = 0;
					for (int l = 0; l < number_of_edges; ++l)
					{
						vertex_index_candidate_0 = from_index_dual[adjacent_vector_indices_h[6*to_index[i] + l]];
						vertex_index_candidate_1 = to_index_dual[adjacent_vector_indices_h[6*to_index[i] + l]];
						retval = in_bool_calculator(vertex_indices, number_of_edges, vertex_index_candidate_0, &check_result);						
						if (check_result == 0)
						{
							vertex_indices[counter] = vertex_index_candidate_0;
							latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
							longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
							++counter;
						}
						retval = in_bool_calculator(vertex_indices, number_of_edges, vertex_index_candidate_1, &check_result);						
						if (check_result == 0)
						{
							vertex_indices[counter] = vertex_index_candidate_1;
							latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
							longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
							++counter;
						}
					}
					if (counter != number_of_edges)
					{
						printf("Problem 15 in TRSK implementation detected.\n");
						exit(1);
					}
					for (int l = 0; l < number_of_edges; ++l)
					{
						if (vertex_indices[l] < 0 || vertex_indices[l] >= NUMBER_OF_DUAL_SCALARS_H)
						{
							printf("Problem 5 in TRSK implementation detected.\n");
							exit(1);
						}
						for (int m = l + 1; m < number_of_edges; ++m)
						{
							if (vertex_indices[l] == vertex_indices[m])
							{
								printf("Problem 23 in TRSK implementation detected.\n");
								exit(1);
							}
						}
					}
					retval = sort_edge_indices(latitude_vertices, longitude_vertices, number_of_edges, indices_resorted);
					for (int l = 0; l < number_of_edges; ++l)
					{
						for (int m = l + 1; m < number_of_edges; ++m)
						{
							if (indices_resorted[l] == indices_resorted[m] || indices_resorted[l] < 0 || indices_resorted[l] > number_of_edges - 1)
							{
								printf("Problem 27 in TRSK implementation detected.\n");
								exit(1);
							}
						}
					}
					for (int l = 0; l < number_of_edges; ++l)
						vertex_indices_resorted[l] = vertex_indices[indices_resorted[l]];
					for (int l = 0; l < number_of_edges; ++l)
					{
						for (int m = 0; m < number_of_edges; ++m)
						{
							if ((from_index_dual[adjacent_vector_indices_h[6*to_index[i] + m]] == vertex_indices_resorted[l] && to_index_dual[adjacent_vector_indices_h[6*to_index[i] + m]] == vertex_indices_resorted[(l + 1)%number_of_edges]) || (to_index_dual[adjacent_vector_indices_h[6*to_index[i] + m]] == vertex_indices_resorted[l] && from_index_dual[adjacent_vector_indices_h[6*to_index[i] + m]] == vertex_indices_resorted[(l + 1)%number_of_edges]))
								edge_indices[l] = adjacent_vector_indices_h[6*to_index[i] + m];
						}
					}
					for (int l = 0; l < number_of_edges; ++l)
					{
						if (edge_indices[l] < 0 || edge_indices[l] >= NUMBER_OF_VECTORS_H)
						{
							printf("Problem 9 in TRSK implementation detected.\n");
							exit(1);
						}
						for (int m = l + 1; m < number_of_edges; ++m)
						{
							if (edge_indices[l] == edge_indices[m])
							{
								printf("Problem 20 in TRSK implementation detected.\n");
								exit(1);
							}
						}
					}
					for (int l = 0; l < number_of_edges; ++l)
					{
						latitude_edges[l] = latitude_vector[edge_indices[l]];
						longitude_edges[l] = longitude_vector[edge_indices[l]];
					}
					check_sum = 0;
					for (int l = 0; l < number_of_edges; ++l)	
					{
						if (l == 0)
							retval = calc_triangle_face(latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[number_of_edges - 1], longitude_edges[number_of_edges - 1], &triangle_0);
						else
							retval = calc_triangle_face(latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[l - 1], longitude_edges[l - 1], &triangle_0);
						retval = calc_triangle_face(latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[l], longitude_edges[l], &triangle_1);
						vector_of_areas[l] = pow(RADIUS + z_scalar[to_index[i]], 2)*(triangle_0 + triangle_1);
						check_sum += vector_of_areas[l];
					}
					if (fabs(check_sum/area[to_index[i]] - 1) > 0.001)
					{
						printf("Problem 32 in TRSK implementation detected.\n");
						exit(1);
					}
					for (int l = 0; l < number_of_edges; ++l)
					{
						if (edge_indices[l] == i)
							last_index = l;
						if (edge_indices[l] == trsk_modified_velocity_indices[10*i + k])
							first_index = (l + 1)%number_of_edges;
							
					}
					if (k == 5 + number_of_edges - 1)
						sum_of_weights = 0;
					else
						double_sum_gen(vector_of_areas, number_of_edges, first_index, last_index, &sum_of_weights);
					if (sum_of_weights < 0 || sum_of_weights/area[from_index[i]] > 1)
					{
						printf("Problem 36 in TRSK implementation detected\n");
						exit(1);
					}
					sum_of_weights = sum_of_weights/area[to_index[i]];
					trsk_modified_weights[10*i + k] = sign_0*(sum_of_weights - 0.5)*sign_1;
				}
			}
			trsk_modified_weights[10*i + k] = -normal_distance_dual[trsk_modified_velocity_indices[10*i + k]]/normal_distance[NUMBER_OF_VECTORS_V + i]*trsk_modified_weights[10*i + k];
		}
		// modification following Gassmann (2018)
		if (to_index[i] < NUMBER_OF_PENTAGONS)
			trsk_modified_curl_indices[10*i + 0] = trsk_modified_velocity_indices[10*i + 8];
		else
			trsk_modified_curl_indices[10*i + 0] = trsk_modified_velocity_indices[10*i + 9];
		trsk_modified_curl_indices[10*i + 1] = trsk_modified_velocity_indices[10*i + 0];
		if (from_index[i] < NUMBER_OF_PENTAGONS)
		{
			trsk_modified_curl_indices[10*i + 2] = trsk_modified_velocity_indices[10*i + 3];
			trsk_modified_curl_indices[10*i + 3] = trsk_modified_velocity_indices[10*i + 5];
			trsk_modified_curl_indices[10*i + 4] = 0;
		}
		else
		{
			trsk_modified_curl_indices[10*i + 2] = trsk_modified_velocity_indices[10*i + 2];
			trsk_modified_curl_indices[10*i + 3] = trsk_modified_velocity_indices[10*i + 4];
			trsk_modified_curl_indices[10*i + 4] = trsk_modified_velocity_indices[10*i + 5];
		}
		if (from_index[i] < NUMBER_OF_PENTAGONS)
			trsk_modified_curl_indices[10*i + 5] = trsk_modified_velocity_indices[10*i + 3];
		else
			trsk_modified_curl_indices[10*i + 5] = trsk_modified_velocity_indices[10*i + 4];
		trsk_modified_curl_indices[10*i + 6] = trsk_modified_velocity_indices[10*i + 5];
		if (to_index[i] < NUMBER_OF_PENTAGONS)
		{
			trsk_modified_curl_indices[10*i + 7] = trsk_modified_velocity_indices[10*i + 8];
			trsk_modified_curl_indices[10*i + 8] = trsk_modified_velocity_indices[10*i + 0];
			trsk_modified_curl_indices[10*i + 9] = 0;
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
		for (int j = 0; j < 10; ++j)
		{
			if (trsk_modified_curl_indices[10*i + j] < 0 || trsk_modified_curl_indices[10*i + j] >= NUMBER_OF_VECTORS_H)
			{
				printf("Problem 30 in TRSK implementation detected.\n");
				exit(1);
			}
		}
        recov_hor_ver_pri_index[4*i + 0] = to_index[i];
        recov_hor_ver_pri_index[4*i + 1] = from_index[i];
        recov_hor_ver_pri_index[4*i + 2] = to_index[i] + NUMBER_OF_VECTORS_PER_LAYER;
        recov_hor_ver_pri_index[4*i + 3] = from_index[i] + NUMBER_OF_VECTORS_PER_LAYER;
    }
    int second_index;
    // doing some more checks
	double value_0, value_1;
	for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			value_0 = normal_distance[NUMBER_OF_VECTORS_V + i]/normal_distance_dual[trsk_modified_velocity_indices[10*i + j]]*trsk_modified_weights[10*i + j];
			if (trsk_modified_velocity_indices[10*i + j] != 0 || (trsk_modified_velocity_indices[10*i + j] == 0 && trsk_modified_weights[10*i + j] != 0))
			{
				second_index = -1;			
				for (int k = 0; k < 10; ++k)
				{
					if (trsk_modified_velocity_indices[10*trsk_modified_velocity_indices[10*i + j] + k] == i && trsk_modified_weights[10*trsk_modified_velocity_indices[10*i + j] + k] != 0)
						second_index = k;
				}
				if (second_index == -1)
				{
					printf("Problem 38 in TRSK implementation detected.\n");
					exit(1);
				}
				value_1 = normal_distance[NUMBER_OF_VECTORS_V + trsk_modified_velocity_indices[10*i + j]]/normal_distance_dual[i]*trsk_modified_weights[10*trsk_modified_velocity_indices[10*i + j] + second_index];
				check_sum = value_0 + value_1;
				if (fabs(check_sum) > 0.001)
				{
					printf("Problem 39 in TRSK implementation detected.\n");
					exit(1);
				}
			}
		}
	}
    free(rel_on_line_dual);
    free(face_of_cell_indices);
	return retval;
}

int set_recov_ver(int recov_ver_index[], int adjacent_vector_indices_h[], double recov_ver_0_pri_weight[], double direction[], double recov_ver_0_curl_weight[], double direction_dual[], double recov_ver_1_pri_weight[], double recov_ver_1_curl_weight[])
{
	int number_of_edges;
	double weight_prefactor;
	for (int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
	{
		number_of_edges = 6;
		weight_prefactor = 2.0/6.0;
		if (i < NUMBER_OF_PENTAGONS)
		{
			number_of_edges = 5;
		    weight_prefactor = 2.0/5.0;
	    }
		for (int j = 0; j < number_of_edges; ++j)
		{
		    recov_ver_index[6*i + j] = adjacent_vector_indices_h[6*i + j];
		    recov_ver_0_pri_weight[6*i + j] = weight_prefactor*cos(direction[recov_ver_index[6*i + j]]);
		    recov_ver_0_curl_weight[6*i + j] = weight_prefactor*cos(direction_dual[recov_ver_index[6*i + j]]);
		    recov_ver_1_pri_weight[6*i + j] = weight_prefactor*sin(direction[recov_ver_index[6*i + j]]);
		    recov_ver_1_curl_weight[6*i + j] = weight_prefactor*sin(direction_dual[recov_ver_index[6*i + j]]);
		}
		if (number_of_edges == 5)
		{
		    recov_ver_index[6*i + 5] = 0;
		    recov_ver_0_pri_weight[6*i + 5] = 0;
		    recov_ver_1_pri_weight[6*i + 5] = 0;
		    recov_ver_0_curl_weight[6*i + 5] = 0;
		    recov_ver_1_curl_weight[6*i + 5] = 0;
		}
	}
	return 0;
}






