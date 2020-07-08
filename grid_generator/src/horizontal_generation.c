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
