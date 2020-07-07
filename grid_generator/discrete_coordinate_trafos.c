/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "discrete_coordinate_trafos.h"
#include "enum_and_typedefs.h"
#include <stdlib.h>
#include <stdio.h>
#include "geos95.h"

int find_angle_change(double angle_0, double angle_1, double *result)
{
    double result_pre = angle_1 - angle_0;
    if (result_pre > M_PI)
        result_pre = result_pre - 2*M_PI;
    if (result_pre < -M_PI)
        result_pre = result_pre + 2*M_PI;
    *result = result_pre;
    return 0;
}

int find_coords_from_triangle_on_face_index(int triangle_on_face_index, int res_id, int *coord_0, int *coord_1, int *coord_0_points_amount)
{
    int check = 1;
    int coord_1_pre = -1;
    int min_index, max_index, points_per_edge;
    int retval = find_points_per_edge(res_id, &points_per_edge);
    max_index = -1;
    while (check == 1)
    {
        ++coord_1_pre;
        *coord_0_points_amount = points_per_edge - coord_1_pre;
        max_index = max_index + *coord_0_points_amount;
        min_index = max_index - (*coord_0_points_amount - 1);
        if (triangle_on_face_index <= max_index && triangle_on_face_index >= min_index)
        {
            *coord_0 = triangle_on_face_index - min_index;
            check = 0;
        }
    }
    *coord_1 = coord_1_pre;
    return retval;
}

int find_triangle_on_face_index_from_coords(int coord_0, int coord_1, int res_id, int *triangle_on_face_index)
{
    int i = 0;
    *triangle_on_face_index = 0;
    int points_per_edge;
    int retval = find_points_per_edge(res_id, &points_per_edge);
    int coord_0_points_amount = points_per_edge;
    while (i < coord_1)
    {
        *triangle_on_face_index += coord_0_points_amount;
        coord_0_points_amount -= 1;
        ++i;
    }
    *triangle_on_face_index += coord_0;
    return retval;
}

int find_triangle_indices_from_h_vector_index(int res_id, int i, int *point_0, int *point_1, int *point_2, int *point_3, int *point_4, int *point_5, int *dual_scalar_on_face_index, int *small_triangle_edge_index, int face_edges[][3], int face_vertices[][3], int edge_vertices [][2], int face_edges_reverse[][3])
{
    int face_index = (i - NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
    int on_face_index = i - (NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
    int triangle_on_face_index = on_face_index/3;
    *small_triangle_edge_index = on_face_index - 3*triangle_on_face_index;
    int retval = find_triangle_edge_points(triangle_on_face_index, face_index, res_id, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, face_vertices, face_edges, face_edges_reverse);
    return retval;
}

int find_triangle_edge_points(int triangle_on_face_index, int face_index, int res_id, int *point_0, int *point_1, int *point_2, int *point_3, int *point_4, int *point_5, int *dual_scalar_on_face_index, int face_vertices[][3], int face_edges[][3], int face_edges_reverse[][3])
{
    int coord_0, coord_1, coord_0_points_amount;
    int retval = find_coords_from_triangle_on_face_index(triangle_on_face_index, res_id, &coord_0, &coord_1, &coord_0_points_amount);
    *dual_scalar_on_face_index = 1 + 2*triangle_on_face_index + coord_1;
    int points_per_edge, scalar_points_per_inner_face;
    retval += find_points_per_edge(res_id, &points_per_edge);
    retval += find_scalar_points_per_inner_face(res_id, &scalar_points_per_inner_face);
    if (coord_1 == 0)
    {
        if (face_edges_reverse[face_index][0] == 0)
            *point_0 = NUMBER_OF_PENTAGONS + face_edges[face_index][0]*points_per_edge + coord_0;
        else
            *point_0 = NUMBER_OF_PENTAGONS + (face_edges[face_index][0] + 1)*points_per_edge - 1 - coord_0;
    }
    else
        *point_0 = NUMBER_OF_PENTAGONS + points_per_edge*NUMBER_OF_EDGES + face_index*scalar_points_per_inner_face + triangle_on_face_index - points_per_edge;
    if (coord_0 == points_per_edge - 1 - coord_1)
    {
        if (face_edges_reverse[face_index][1] == 0)
            *point_1 = NUMBER_OF_PENTAGONS + face_edges[face_index][1]*points_per_edge + coord_1;
        else
            *point_1 = NUMBER_OF_PENTAGONS + (face_edges[face_index][1] + 1)*points_per_edge - 1 - coord_1;
    }
    else
        *point_1 = NUMBER_OF_PENTAGONS + points_per_edge*NUMBER_OF_EDGES + face_index*scalar_points_per_inner_face + triangle_on_face_index - coord_1;
    if (coord_0 == 0)
    {
        if (face_edges_reverse[face_index][2] == 0)
            *point_2 = NUMBER_OF_PENTAGONS + (face_edges[face_index][2] + 1)*points_per_edge - 1 - coord_1;
        else
            *point_2 = NUMBER_OF_PENTAGONS + face_edges[face_index][2]*points_per_edge + coord_1;
    }
    else
        *point_2 = NUMBER_OF_PENTAGONS + points_per_edge*NUMBER_OF_EDGES + face_index*scalar_points_per_inner_face + triangle_on_face_index - 1 - coord_1;
    if (coord_1 == 0)
    {
        if (coord_0 == 0)
            *point_3 = face_vertices[face_index][0];
        else
        {
            if (face_edges_reverse[face_index][0] == 0)
                *point_3 = *point_0 - 1;
            else
                *point_3 = *point_0 + 1;
        }
    }
    else if (coord_0 == 0)
    {
        if (face_edges_reverse[face_index][2] == 0)
            *point_3 = *point_2 + 1;
        else
            *point_3 = *point_2 - 1;
    }
    else
        *point_3 = *point_0 - 1;
    *point_4 = -1;
    *point_5 = -1;
    if (coord_0 == coord_0_points_amount - 1)
    {
        if (coord_1 == 0)
            *point_4 = face_vertices[face_index][1];
        else
        {
            if (face_edges_reverse[face_index][1] == 0)
                *point_4 = *point_1 - 1;
            else
                *point_4 = *point_1 + 1;
        }
        if (coord_1 == points_per_edge - 1)
            *point_5 = face_vertices[face_index][2];
    }
    return 0;
}

int find_triangle_on_face_index_from_dual_scalar_on_face_index(int dual_scalar_on_face_index, int res_id, int *triangle_on_face_index, int *points_downwards, int *special_case_bool, int *last_triangle_bool)
{
    int value_found = 0;
    int retval = 0;
    int triangle_on_face_index_pre, coord_0_pre, coord_1_pre, coord_0_points_amount_pre, dual_scalar_on_face_index_0, dual_scalar_on_face_index_1, dual_scalar_on_face_index_2, dual_scalar_on_face_index_3, points_per_edge;
    triangle_on_face_index_pre = -1;
    retval += find_points_per_edge(res_id, &points_per_edge);
    while (value_found == 0)
    {
        dual_scalar_on_face_index_2 = -1;
        dual_scalar_on_face_index_3 = -1;
        ++triangle_on_face_index_pre;
        retval += find_coords_from_triangle_on_face_index(triangle_on_face_index_pre, res_id, &coord_0_pre, &coord_1_pre, &coord_0_points_amount_pre);
        dual_scalar_on_face_index_0 = 2*triangle_on_face_index_pre + 1 + coord_1_pre;
        dual_scalar_on_face_index_1 = dual_scalar_on_face_index_0 - 1;
        if (coord_0_pre == coord_0_points_amount_pre - 1)
        {
            dual_scalar_on_face_index_2 = dual_scalar_on_face_index_0 + 1;
            if (coord_1_pre == points_per_edge - 1)
                dual_scalar_on_face_index_3 = dual_scalar_on_face_index_2 + 1;
        }
        if (dual_scalar_on_face_index == dual_scalar_on_face_index_0)
        {
            *points_downwards = 1;
            *special_case_bool = 0;
            *last_triangle_bool = 0;
            value_found = 1;
        }
        if (dual_scalar_on_face_index == dual_scalar_on_face_index_1)
        {
            *points_downwards = 0;
            *special_case_bool = 0;
            *last_triangle_bool = 0;
            value_found = 1;
        }
        if (dual_scalar_on_face_index == dual_scalar_on_face_index_2)
        {
            *points_downwards = 0;
            *special_case_bool = 1;
            *last_triangle_bool = 0;
            value_found = 1;
        }
        if (dual_scalar_on_face_index == dual_scalar_on_face_index_3)
        {
            *points_downwards = 0;
            *special_case_bool = 0;
            *last_triangle_bool = 1;
            value_found = 1;
        }
    }
    *triangle_on_face_index = triangle_on_face_index_pre;
    return retval;
}

int find_triangle_edge_points_from_dual_scalar_on_face_index(int dual_scalar_on_face_index, int face_index, int res_id, int *point_0, int *point_1, int *point_2, int face_vertices[][3], int face_edges[][3], int face_edges_reverse[][3])
{
    int points_downwards, special_case_bool, last_triangle_bool;
    int triangle_on_face_index, rhombuspoint_0, rhombuspoint_1, rhombuspoint_2, rhombuspoint_3, coord_0, coord_1, coord_0_points_amount, points_per_edge, dump, addpoint_0, addpoint_1;
    int retval = find_triangle_on_face_index_from_dual_scalar_on_face_index(dual_scalar_on_face_index, res_id, &triangle_on_face_index, &points_downwards, &special_case_bool, &last_triangle_bool);
    retval += find_coords_from_triangle_on_face_index(triangle_on_face_index, res_id, &coord_0, &coord_1, &coord_0_points_amount);
    retval += find_points_per_edge(res_id, &points_per_edge);
    retval += find_triangle_edge_points(triangle_on_face_index, face_index, res_id, &rhombuspoint_0, &rhombuspoint_1, &rhombuspoint_2, &rhombuspoint_3, &addpoint_0, &addpoint_1, &dump, face_vertices, face_edges, face_edges_reverse);
    if (points_downwards == 1)
    {
        *point_0 = rhombuspoint_0;
        *point_1 = rhombuspoint_1;
        *point_2 = rhombuspoint_2;
    }
    else
    {
        if (coord_0 == coord_0_points_amount - 1)
        {
            if (coord_1 == points_per_edge - 1)
            {
                if (last_triangle_bool == 1)
                {
                    *point_0 = rhombuspoint_2;
                    *point_1 = rhombuspoint_1;
                    *point_2 = addpoint_1;   
                }
                else if (special_case_bool == 1)
                {
                    *point_0 = rhombuspoint_0;
                    *point_1 = addpoint_0;
                    *point_2 = rhombuspoint_1;   
                }
                else
                {
                    *point_0 = rhombuspoint_3;
                    *point_1 = rhombuspoint_0;
                    *point_2 = rhombuspoint_2;  
                }
            }
            else
            {
                if (special_case_bool == 1)
                {
                    *point_0 = rhombuspoint_0;
                    *point_1 = addpoint_0;
                    *point_2 = rhombuspoint_1;  
                }
                else
                {
                    *point_0 = rhombuspoint_3;
                    *point_1 = rhombuspoint_0;
                    *point_2 = rhombuspoint_2;  
                }
            }
        }
        else
        {
            *point_0 = rhombuspoint_3;
            *point_1 = rhombuspoint_0;
            *point_2 = rhombuspoint_2;
        }
    }
    return retval;
}

int find_points_per_edge(int res_id, int *points_per_edge)
{
    *points_per_edge = (int) (pow(2, res_id) - 1);
    return 0;
}

int find_scalar_points_per_inner_face(int res_id, int *scalar_points_per_inner_face)
{
    *scalar_points_per_inner_face = (int) (0.5*(pow(2, res_id) - 2)*(pow(2, res_id) - 1));
    return 0;
}

int upscale_scalar_point(int res_id, int old_index, int *new_index)
{
    int edge_index, face_index;
    int points_per_edge, on_edge_index, scalar_points_per_inner_face, on_face_index, coord_0, coord_1, coord_0_points_amount;
    int retval = find_points_per_edge(res_id, &points_per_edge);
    retval += find_scalar_points_per_inner_face(res_id, &scalar_points_per_inner_face);
    if (old_index < NUMBER_OF_PENTAGONS)
        *new_index = old_index;
    else if (old_index < NUMBER_OF_PENTAGONS + NUMBER_OF_EDGES*points_per_edge)
    {
        edge_index = (old_index - NUMBER_OF_PENTAGONS)/points_per_edge;
        on_edge_index = old_index - (NUMBER_OF_PENTAGONS + edge_index*points_per_edge);
        *new_index = NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE + pow(2, RES_ID - res_id)*(on_edge_index + 1) - 1;
    }
    else
    {
        face_index = (old_index - (NUMBER_OF_PENTAGONS + NUMBER_OF_EDGES*points_per_edge))/scalar_points_per_inner_face;
        on_face_index = old_index - (NUMBER_OF_PENTAGONS + NUMBER_OF_EDGES*points_per_edge + face_index*scalar_points_per_inner_face);
        retval = find_coords_from_triangle_on_face_index(on_face_index + points_per_edge, res_id, &coord_0, &coord_1, &coord_0_points_amount);
        coord_0 = (coord_0 + 1)*pow(2, RES_ID - res_id) - 1;
        coord_1 = coord_1*pow(2, RES_ID - res_id);
        retval += find_triangle_on_face_index_from_coords(coord_0, coord_1, RES_ID, &on_face_index);
        *new_index = NUMBER_OF_PENTAGONS + NUMBER_OF_EDGES*POINTS_PER_EDGE + face_index*SCALAR_POINTS_PER_INNER_FACE + on_face_index - POINTS_PER_EDGE;
    }
    return retval;
}

int write_scalar_coordinates(int edgepoint_0, int edgepoint_1, int edgepoint_2, int point_0, int point_1, int point_2, int points_upwards, double x_unity[], double y_unity[], double z_unity[], double latitude_scalar[], double longitude_scalar[])
{
    double x_res, y_res, z_res, lat_res, lon_res;
    int retval = find_between_point(x_unity[edgepoint_0], y_unity[edgepoint_0], z_unity[edgepoint_0], x_unity[edgepoint_1], y_unity[edgepoint_1], z_unity[edgepoint_1], 0.5, &x_res, &y_res, &z_res);
    retval += normalize_cartesian(x_res, y_res, z_res, &x_res, &y_res, &z_res);
    if (points_upwards == 1)
    {
        x_unity[point_0] = x_res;
        y_unity[point_0] = y_res;
        z_unity[point_0] = z_res;
    }
    else
    {
        x_unity[point_1] = x_res;
        y_unity[point_1] = y_res;
        z_unity[point_1] = z_res;
    }
    retval = find_geos(x_res, y_res, z_res, &lat_res, &lon_res);
    if (points_upwards == 1)
    {
        latitude_scalar[point_0] = lat_res;
        longitude_scalar[point_0] = lon_res;
    }
    else
    {
        latitude_scalar[point_1] = lat_res;
        longitude_scalar[point_1] = lon_res;
    }
    retval += find_between_point(x_unity[edgepoint_1], y_unity[edgepoint_1], z_unity[edgepoint_1], x_unity[edgepoint_2], y_unity[edgepoint_2], z_unity[edgepoint_2], 0.5, &x_res, &y_res, &z_res);
    retval += normalize_cartesian(x_res, y_res, z_res, &x_res, &y_res, &z_res);
    if (points_upwards == 1)
    {
        x_unity[point_1] = x_res;
        y_unity[point_1] = y_res;
        z_unity[point_1] = z_res;
    }
    else
    {
        x_unity[point_2] = x_res;
        y_unity[point_2] = y_res;
        z_unity[point_2] = z_res;
    }
    retval += find_geos(x_res, y_res, z_res, &lat_res, &lon_res);
    if (points_upwards == 1)
    {
        latitude_scalar[point_1] = lat_res;
        longitude_scalar[point_1] = lon_res;
    }
    else
    {
        latitude_scalar[point_2] = lat_res;
        longitude_scalar[point_2] = lon_res;
    }
    retval += find_between_point(x_unity[edgepoint_2], y_unity[edgepoint_2], z_unity[edgepoint_2], x_unity[edgepoint_0], y_unity[edgepoint_0], z_unity[edgepoint_0], 0.5, &x_res, &y_res, &z_res);
    retval += normalize_cartesian(x_res, y_res, z_res, &x_res, &y_res, &z_res);
    if (points_upwards == 1)
    {
        x_unity[point_2] = x_res;
        y_unity[point_2] = y_res;
        z_unity[point_2] = z_res;
    }
    else
    {
        x_unity[point_0] = x_res;
        y_unity[point_0] = y_res;
        z_unity[point_0] = z_res;
    }
    retval += find_geos(x_res, y_res, z_res, &lat_res, &lon_res);
    if (points_upwards == 1)
    {
        latitude_scalar[point_2] = lat_res;
        longitude_scalar[point_2] = lon_res;
    }
    else
    {
        latitude_scalar[point_0] = lat_res;
        longitude_scalar[point_0] = lon_res;
    }
    return retval;
}

int find_triangles_per_face(int res_id, int *number_of_triangles_per_face)
{
    *number_of_triangles_per_face = (int) (pow(4, res_id));
    return 0;
}

int find_v_vector_indices_for_dual_scalar_z(int from_index[], int to_index[], int vorticity_indices_pre[], int dual_scalar_h_index, int index_vector_for_dual_scalar_z[])
{
	int counter = 0;
	int check_result, retval;
	index_vector_for_dual_scalar_z[0] = -1;
	index_vector_for_dual_scalar_z[1] = -1;
	index_vector_for_dual_scalar_z[2] = -1;
	for (int k = 0; k < 3; ++k)
	{
		retval = in_bool_calculator(index_vector_for_dual_scalar_z, 3, from_index[vorticity_indices_pre[3*dual_scalar_h_index + k]], &check_result);
		if (check_result == 0)
		{
			index_vector_for_dual_scalar_z[counter] = from_index[vorticity_indices_pre[3*dual_scalar_h_index + k]];
			counter++;
		}
		retval = in_bool_calculator(index_vector_for_dual_scalar_z, 3, to_index[vorticity_indices_pre[3*dual_scalar_h_index + k]], &check_result);
		if (check_result == 0)
		{
			index_vector_for_dual_scalar_z[counter] = to_index[vorticity_indices_pre[3*dual_scalar_h_index + k]];
			counter++;
		}
	}
	if (counter != 3)
	{
		printf("Error in function find_v_vector_indices_for_dual_scalar_z.\n");
		exit(1);
	}
	if (retval != 0)
		return 1;
	else
		return 0;
}














