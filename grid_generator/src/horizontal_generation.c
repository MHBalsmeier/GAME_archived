#include "enum.h"
#include "grid_generator.h"
#include <stdlib.h>
#include <stdio.h>
#include "geos95.h"

int generate_horizontal_generators(double latitude_ico[], double longitude_ico[], double latitude_scalar[], double longitude_scalar[], double x_unity[], double y_unity[], double z_unity[], int face_edges_reverse[][3], int face_edges[][3], int face_vertices[][3])
{
	int base_index_down_triangles, base_index_old, test_index, last_triangle_bool, old_triangle_on_line_index, base_index_up_triangles, points_downwards, points_upwards, dump, points_per_edge, edgepoint_0, edgepoint_1, edgepoint_2, number_of_triangles_per_face, retval, point_0, point_1, point_2, dual_scalar_on_face_index, coord_0, coord_1, triangle_on_face_index, coord_0_points_amount;
	double x_res, y_res, z_res;
	for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
	{
	    upscale_scalar_point(RES_ID, i, &test_index);
	    if (test_index != i)
	        printf("problem with upscale_scalar_point detected\n");
	}
	for (int i = 0; i < NUMBER_OF_PENTAGONS; ++i)
	{
	    latitude_scalar[i] = latitude_ico[i];
	    longitude_scalar[i] = longitude_ico[i];
	    retval = find_global_normal(latitude_ico[i], longitude_ico[i], &x_res, &y_res, &z_res);
	    x_unity[i] = x_res;
	    y_unity[i] = y_res;
	    z_unity[i] = z_res;
	}
	for (int i = 0; i < NUMBER_OF_BASIC_TRIANGLES; ++i)
	{
	    for (int j = 0; j < RES_ID; ++j)
	    {
	        retval = find_triangles_per_face(j, &number_of_triangles_per_face);
	        for (int k = 0; k < number_of_triangles_per_face; ++k)
	        {
	            if (j == 0)
	            {
	                dual_scalar_on_face_index = 1;
	                retval = find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index, i, j + 1, &point_0, &point_1, &point_2, face_vertices, face_edges, face_edges_reverse);
	                retval += upscale_scalar_point(j + 1, point_0, &point_0);
	                retval += upscale_scalar_point(j + 1, point_1, &point_1);
	                retval += upscale_scalar_point(j + 1, point_2, &point_2);
	                points_upwards = 1;
	                retval = write_scalar_coordinates(face_vertices[i][0], face_vertices[i][1], face_vertices[i][2], point_0, point_1, point_2, points_upwards, x_unity, y_unity, z_unity, latitude_scalar, longitude_scalar);
	            }
	            else
	            {
	                retval = find_triangle_edge_points_from_dual_scalar_on_face_index(k, i, j, &edgepoint_0, &edgepoint_1, &edgepoint_2, face_vertices, face_edges, face_edges_reverse);
	                retval += find_triangle_on_face_index_from_dual_scalar_on_face_index(k, j, &triangle_on_face_index, &points_downwards, &dump, &last_triangle_bool);
	                retval += find_coords_from_triangle_on_face_index(triangle_on_face_index, j, &coord_0, &coord_1, &coord_0_points_amount);
	                retval += find_points_per_edge(j, &points_per_edge);
	                base_index_old = 0;
	                base_index_down_triangles = 0;
	                base_index_up_triangles = base_index_down_triangles + 4*points_per_edge + 3;
	                for (int l = 0; l < coord_1; ++l)
	                {
	                    coord_0_points_amount = points_per_edge - l;
	                    base_index_old += 2*coord_0_points_amount + 1;
	                    base_index_down_triangles += 4*(2*coord_0_points_amount + 1);
	                    base_index_up_triangles = base_index_down_triangles + 4*(points_per_edge - l) + 3;
	                }
	                if (last_triangle_bool == 1)
	                {
	                    base_index_old += 3;
	                    base_index_down_triangles += 12;
	                    base_index_up_triangles = base_index_down_triangles + 3;
	                }
	                old_triangle_on_line_index = k - base_index_old;
	                if (points_downwards == 0)
	                    dual_scalar_on_face_index = base_index_down_triangles + 1 + 2*old_triangle_on_line_index;
	                else
	                    dual_scalar_on_face_index = base_index_up_triangles + 2*old_triangle_on_line_index;
	                retval = find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index, i, j + 1, &point_0, &point_1, &point_2, face_vertices, face_edges, face_edges_reverse);
	                retval += upscale_scalar_point(j, edgepoint_0, &edgepoint_0);
	                retval += upscale_scalar_point(j, edgepoint_1, &edgepoint_1);
	                retval += upscale_scalar_point(j, edgepoint_2, &edgepoint_2);
	                retval += upscale_scalar_point(j + 1, point_0, &point_0);
	                retval += upscale_scalar_point(j + 1, point_1, &point_1);
	                retval += upscale_scalar_point(j + 1, point_2, &point_2);
	                points_upwards = 1;
	                if (points_downwards == 1)
	                    points_upwards = 0;
	                retval = write_scalar_coordinates(edgepoint_0, edgepoint_1, edgepoint_2, point_0, point_1, point_2, points_upwards, x_unity, y_unity, z_unity, latitude_scalar, longitude_scalar);
	            }
	        }
	    }
	}
	free(latitude_ico);
	free(longitude_ico);
	return 0;
}

int calc_cell_face_unity(double pent_hex_face_unity_sphere[], double latitude_scalar_dual[], double longitude_scalar_dual[], int adjacent_vector_indices_h [], int vorticity_indices_pre [])
{
    int check_0, check_1, check_2, retval, counter;
    for (int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
    {
        if (i < NUMBER_OF_PENTAGONS)
        {
            double *lat_points = malloc(5*sizeof(double));
            double *lon_points = malloc(5*sizeof(double));
            int *cell_vector_indices = malloc(5*sizeof(int));
            for (int j = 0; j < 5; ++j)
                cell_vector_indices[j] = adjacent_vector_indices_h[6*i + j];
            counter = 0;
            for (int j = 0; j < NUMBER_OF_DUAL_VECTORS_V; ++j)
            {
                retval = in_bool_calculator(cell_vector_indices, 5, vorticity_indices_pre[3*j + 0], &check_0);
                retval = in_bool_calculator(cell_vector_indices, 5, vorticity_indices_pre[3*j + 1], &check_1);
                retval = in_bool_calculator(cell_vector_indices, 5, vorticity_indices_pre[3*j + 2], &check_2);
                if (check_0 == 1 || check_1 == 1 || check_2 == 1)
                {
                    lat_points[counter] = latitude_scalar_dual[j];
                    lon_points[counter] = longitude_scalar_dual[j];
                    counter++;
                }
            }
            retval = calc_spherical_polygon_face(lat_points, lon_points, 5, &pent_hex_face_unity_sphere[i]);
            free(lat_points);
            free(lon_points);
            free(cell_vector_indices);
        }
        else
        {
            double *lat_points = malloc(6*sizeof(double));
            double *lon_points = malloc(6*sizeof(double));
            int *cell_vector_indices = malloc(6*sizeof(int));
            for (int j = 0; j < 6; ++j)
                cell_vector_indices[j] = adjacent_vector_indices_h[6*i + j];
            counter = 0;
            for (int j = 0; j < NUMBER_OF_DUAL_VECTORS_V; ++j)
            {
                retval = in_bool_calculator(cell_vector_indices, 6, vorticity_indices_pre[3*j + 0], &check_0);
                retval = in_bool_calculator(cell_vector_indices, 6, vorticity_indices_pre[3*j + 1], &check_1);
                retval = in_bool_calculator(cell_vector_indices, 6, vorticity_indices_pre[3*j + 2], &check_2);
                if (check_0 == 1 || check_1 == 1 || check_2 == 1)
                {
                    lat_points[counter] = latitude_scalar_dual[j];
                    lon_points[counter] = longitude_scalar_dual[j];
                    counter++;
                }
            }
            retval = calc_spherical_polygon_face(lat_points, lon_points, 6, &pent_hex_face_unity_sphere[i]);
            free(lat_points);
            free(lon_points);
            free(cell_vector_indices);
        }
    }
    double pent_hex_sum_unity_sphere = 0;
    double pent_hex_avg_unity_sphere_ideal = 4*M_PI/NUMBER_OF_SCALARS_H;
    for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
    {
        pent_hex_sum_unity_sphere += pent_hex_face_unity_sphere[i];
        if (pent_hex_face_unity_sphere[i] <= 0)
		{
            printf("pent_hex_face_unity_sphere contains a non-positive value.\n");
			exit(1);
		}
        if (fabs(pent_hex_face_unity_sphere[i]/pent_hex_avg_unity_sphere_ideal - 1) > 0.4)
		{
            printf("Pentagons and hexagons on unity sphere have significantly different surfaces.\n");
			exit(1);
		}
    }
    if (fabs(pent_hex_sum_unity_sphere/(4*M_PI) - 1) > 1e-12)
	{
        printf("Sum of faces of pentagons and hexagons on unity sphere does not match face of unit sphere.\n");
		exit(1);
	}
	return retval;
}

int calc_triangle_face_unity(double triangle_face_unit_sphere[], double latitude_scalar[], double longitude_scalar[], int face_edges[][3], int face_edges_reverse[][3], int face_vertices[][3], int edge_vertices[][2])
{
	int retval, dual_scalar_index, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, small_triangle_edge_index, coord_0_points_amount, coord_0, coord_1, face_index, on_face_index, triangle_on_face_index;
	double triangle_face;
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        if (i >= NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))
        {
            retval = find_triangle_indices_from_h_vector_index(RES_ID, i, &point_0, &point_1, &point_2, &point_3, &point_4, &point_5, &dual_scalar_on_face_index, &small_triangle_edge_index, face_edges, face_vertices, edge_vertices, face_edges_reverse);
            face_index = (i - NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
            on_face_index = i - (NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
            triangle_on_face_index = on_face_index/3;
            retval += find_coords_from_triangle_on_face_index(triangle_on_face_index, RES_ID, &coord_0, &coord_1, &coord_0_points_amount);
            dual_scalar_index = dual_scalar_on_face_index + face_index*NUMBER_OF_TRIANGLES/NUMBER_OF_BASIC_TRIANGLES;
            retval += calc_triangle_face(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_2], longitude_scalar[point_2], &triangle_face);
            triangle_face_unit_sphere[dual_scalar_index] = triangle_face;
            retval += calc_triangle_face(latitude_scalar[point_3], longitude_scalar[point_3], latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_2], longitude_scalar[point_2], &triangle_face);
            triangle_face_unit_sphere[dual_scalar_index - 1] = triangle_face;
            if (coord_0 == coord_0_points_amount - 1)
            {
                retval += calc_triangle_face(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_4], longitude_scalar[point_4], latitude_scalar[point_1], longitude_scalar[point_1], &triangle_face);
                triangle_face_unit_sphere[dual_scalar_index + 1] = triangle_face;
                if (coord_1 == POINTS_PER_EDGE - 1)
                {
                    retval += calc_triangle_face(latitude_scalar[point_2], longitude_scalar[point_2], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_5], longitude_scalar[point_5], &triangle_face);
                    triangle_face_unit_sphere[dual_scalar_index + 2] = triangle_face;
                }
            }
        }
    }
    double triangle_sum_unit_sphere = 0;
    double triangle_avg_unit_sphere_ideal = 4*M_PI/NUMBER_OF_TRIANGLES;
    for (int i = 0; i < NUMBER_OF_DUAL_SCALARS_H; ++i)
    {
        triangle_sum_unit_sphere += triangle_face_unit_sphere[i];
        if (triangle_face_unit_sphere[i] <= 0)
		{
            printf("triangle_face_unit_sphere contains a non-positive value.\n");
			exit(1);
		}
        if (fabs(triangle_face_unit_sphere[i]/triangle_avg_unit_sphere_ideal - 1) > 0.4)
		{
            printf("Triangles on unit sphere have significantly different surfaces.\n");
			exit(1);
		}
    }
    if (fabs(triangle_sum_unit_sphere/(4*M_PI) - 1) > 1e-12)
	{
        printf("Sum of faces of triangles on unit sphere does not match face of unit sphere.\n");
		exit(1);
	}
    return 0;
}


int set_vector_h_doubles(int from_index[], int to_index[], double latitude_scalar[], double longitude_scalar[], double latitude_vector[], double longitude_vector[], double direction[])
{
	int retval;
	double x_point_0, y_point_0, z_point_0, x_point_1, y_point_1, z_point_1, x_res, y_res, z_res, lat_res, lon_res;
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        retval = find_global_normal(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], &x_point_0, &y_point_0, &z_point_0);
        retval = find_global_normal(latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], &x_point_1, &y_point_1, &z_point_1);
        retval = find_between_point(x_point_0, y_point_0, z_point_0, x_point_1, y_point_1, z_point_1, 0.5, &x_res, &y_res, &z_res);
        retval = find_geos(x_res, y_res, z_res, &lat_res, &lon_res);
        latitude_vector[i] = lat_res;
        longitude_vector[i] = lon_res;
        direction[i] = find_geodetic_direction(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], 0.5);
    }
	return retval;
}

int set_from_to_index(int from_index[], int to_index[], int face_edges[][3], int face_edges_reverse[][3], int face_vertices[][3], int edge_vertices[][2])
{
	int edge_index, on_edge_index, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, small_triangle_edge_index, retval;
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        if (i < NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))
        {
            edge_index = i/(POINTS_PER_EDGE + 1);
            on_edge_index = i - edge_index*(POINTS_PER_EDGE + 1);
            if(on_edge_index == 0)
            {
                from_index[i] = edge_vertices[edge_index][0];
                to_index[i] = NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE;
            }
            else if (on_edge_index == POINTS_PER_EDGE)
            {
                from_index[i] = NUMBER_OF_PENTAGONS + (edge_index + 1)*POINTS_PER_EDGE - 1;
                to_index[i] = edge_vertices[edge_index][1];
            }
            else
            {
                from_index[i] = NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE + on_edge_index - 1;
                to_index[i] = NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE + on_edge_index;
            }
        }
        else
        {
            retval = find_triangle_indices_from_h_vector_index(RES_ID, i, &point_0, &point_1, &point_2, &point_3, &point_4, &point_5, &dual_scalar_on_face_index, &small_triangle_edge_index, face_edges, face_vertices, edge_vertices, face_edges_reverse);
            if (small_triangle_edge_index == 0)
            {
                from_index[i] = point_0;
                to_index[i] = point_2;
            }
            if (small_triangle_edge_index == 1)
            {
                from_index[i] = point_0;
                to_index[i] = point_1;
            }
            if (small_triangle_edge_index == 2)
            {
                from_index[i] = point_2;
                to_index[i] = point_1;
            }
        }
    }
    return retval;
}

int set_scalar_h_dual_coords(double latitude_scalar_dual[], double longitude_scalar_dual[], double latitude_scalar[], double longitude_scalar[], int face_edges[][3], int face_edges_reverse[][3], int face_vertices[][3], int edge_vertices[][2])
{
	double lat_res, lon_res;
	int retval, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, small_triangle_edge_index, dual_scalar_index, coord_0, coord_1, coord_0_points_amount, face_index, on_face_index, triangle_on_face_index;
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        if (i >= NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))
        {
            retval = find_triangle_indices_from_h_vector_index(RES_ID, i, &point_0, &point_1, &point_2, &point_3, &point_4, &point_5, &dual_scalar_on_face_index, &small_triangle_edge_index, face_edges, face_vertices, edge_vertices, face_edges_reverse);
            face_index = (i - NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
            on_face_index = i - (NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
            triangle_on_face_index = on_face_index/3;
            retval += find_coords_from_triangle_on_face_index(triangle_on_face_index, RES_ID, &coord_0, &coord_1, &coord_0_points_amount);
            dual_scalar_index = dual_scalar_on_face_index + face_index*NUMBER_OF_TRIANGLES/NUMBER_OF_BASIC_TRIANGLES;
            retval += find_voronoi_center_sphere(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_2], longitude_scalar[point_2], &lat_res, &lon_res);
            latitude_scalar_dual[dual_scalar_index] = lat_res;
            longitude_scalar_dual[dual_scalar_index] = lon_res;
	        retval += find_voronoi_center_sphere(latitude_scalar[point_3], longitude_scalar[point_3], latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_2], longitude_scalar[point_2], &lat_res, &lon_res);
            latitude_scalar_dual[dual_scalar_index - 1] = lat_res;
            longitude_scalar_dual[dual_scalar_index - 1] = lon_res;
            if (coord_0 == coord_0_points_amount - 1)
            {
                retval += find_voronoi_center_sphere(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_4], longitude_scalar[point_4], latitude_scalar[point_1], longitude_scalar[point_1], &lat_res, &lon_res);
                latitude_scalar_dual[dual_scalar_index + 1] = lat_res;
                longitude_scalar_dual[dual_scalar_index + 1] = lon_res;
                if (coord_1 == POINTS_PER_EDGE - 1)
                {
                    retval += find_voronoi_center_sphere(latitude_scalar[point_2], longitude_scalar[point_2], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_5], longitude_scalar[point_5], &lat_res, &lon_res);
                    latitude_scalar_dual[dual_scalar_index + 2] = lat_res;
                    longitude_scalar_dual[dual_scalar_index + 2] = lon_res;
                }
            }
        }
    }
	return 0;
}












