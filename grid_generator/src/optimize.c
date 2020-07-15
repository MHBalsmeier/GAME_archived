/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
The Lloyd algorithm is implemented here.
*/

#include "enum.h"
#include "grid_generator.h"
#include <stdlib.h>
#include <stdio.h>
#include "geos95.h"

int find_cell_cgs(double [], double [], double [], double [], int [], int [], int []);

int optimize_to_scvt(double latitude_scalar[], double longitude_scalar[], double latitude_scalar_dual[], double longitude_scalar_dual[], int n_iterations, int face_edges[][3], int face_edges_reverse[][3], int face_vertices[][3], int edge_vertices[][2], int adjacent_vector_indices_h[], int from_index_dual[], int to_index_dual[])
{
	int retval = 0;
	for (int i = 0; i < n_iterations; ++i)
	{
    	retval = set_scalar_h_dual_coords(latitude_scalar_dual, longitude_scalar_dual, latitude_scalar, longitude_scalar, face_edges, face_edges_reverse, face_vertices, edge_vertices);
    	retval = find_cell_cgs(latitude_scalar, longitude_scalar, latitude_scalar_dual, longitude_scalar_dual, adjacent_vector_indices_h, from_index_dual, to_index_dual);
    	printf("Optimizing grid - iteration %d completed.\n", i + 1);
	}
	return retval;
}

int find_cell_cgs(double latitude_scalar[], double longitude_scalar[], double latitude_scalar_dual[], double longitude_scalar_dual[], int adjacent_vector_indices_h[], int from_index_dual[], int to_index_dual[])
{
	int number_of_edges, counter, vertex_index_candidate_0, vertex_index_candidate_1, retval, check_result;
	double lat_res, lon_res, x_res, y_res, z_res, triangle_unity_face, x_0, y_0, z_0, x_1, y_1, z_1, x_2, y_2, z_2, lat_0, lon_0, lat_1, lon_1, lat_2, lon_2;
	for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
	{
		number_of_edges = 6;
		if (i < NUMBER_OF_PENTAGONS)
			number_of_edges = 5;
		double latitude_vertices[number_of_edges];
		double longitude_vertices[number_of_edges];
		int vertex_indices[number_of_edges];
		int vertex_indices_resorted[number_of_edges];
		int indices_resorted[number_of_edges];
		for (int l = 0; l < number_of_edges; ++l)
			vertex_indices[l] = -1;
		counter = 0;
		for (int j = 0; j < number_of_edges; ++j)
		{
			vertex_index_candidate_0 = from_index_dual[adjacent_vector_indices_h[6*i + j]];
			vertex_index_candidate_1 = to_index_dual[adjacent_vector_indices_h[6*i + j]];
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
			printf("Trouble in find_cell_cgs detected.\n");
		sort_edge_indices(latitude_vertices, longitude_vertices, number_of_edges, indices_resorted);
		for (int j = 0; j < number_of_edges; ++j)
			vertex_indices_resorted[j] = vertex_indices[indices_resorted[j]];
		x_res = 0;
		y_res = 0;
		z_res = 0;
		for (int j = 0; j < number_of_edges; ++j)
		{
			lat_0 = latitude_scalar[i];
			lon_0 = longitude_scalar[i];
			lat_1 = latitude_scalar_dual[vertex_indices_resorted[j]];
			lon_1 = longitude_scalar_dual[vertex_indices_resorted[j]];
			lat_2 = latitude_scalar_dual[vertex_indices_resorted[(j + 1)%number_of_edges]];
			lon_2 = longitude_scalar_dual[vertex_indices_resorted[(j + 1)%number_of_edges]];
			find_global_normal(lat_0, lon_0, &x_0, &y_0, &z_0);
			find_global_normal(lat_1, lon_1, &x_1, &y_1, &z_1);
			find_global_normal(lat_2, lon_2, &x_2, &y_2, &z_2);
			calc_triangle_face(lat_0, lon_0, lat_1, lon_1, lat_2, lon_2, &triangle_unity_face);
			x_res += triangle_unity_face*(x_0 + x_1 + x_2);
			y_res += triangle_unity_face*(y_0 + y_1 + y_2);
			z_res += triangle_unity_face*(z_0 + z_1 + z_2);
		}
		find_geos(x_res, y_res, z_res, &lat_res, &lon_res);
		latitude_scalar[i] = lat_res;
		longitude_scalar[i] = lon_res;
	}
	return retval;
}





