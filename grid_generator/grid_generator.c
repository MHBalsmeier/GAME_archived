/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "discrete_coordinate_trafos.h"
#include "enum.h"
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <math.h>
#include "geos95.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define OMEGA (7.292115e-5)
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define C_P 1005.0
#define P_0 100000.0
#define GRAVITY_MEAN_SFC_ABS 9.80616

const double TOA = 30000.0;
// code fails for ORO_ID = 2 (steep orography test)
const double ORTH_CRITERION_DEG = 89.99;

int main(int argc, char *argv[])
{
    int ORO_ID;
   	ORO_ID = strtod(argv[1], NULL);
	if (NUMBER_OF_ORO_LAYERS >= NUMBER_OF_LAYERS)
	{
		printf("It is NUMBER_OF_ORO_LAYERS >= NUMBER_OF_LAYERS.\n");
		exit(1);
	}
    int OUTPUT_FILE_LENGTH = 100;
    char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE_PRE, "nc_files/B%dL%dT%d_O%d_OL%d.nc", RES_ID, NUMBER_OF_LAYERS, (int) TOA, ORO_ID, NUMBER_OF_ORO_LAYERS);
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "nc_files/B%dL%dT%d_O%d_OL%d.nc", RES_ID, NUMBER_OF_LAYERS, (int) TOA, ORO_ID, NUMBER_OF_ORO_LAYERS);
    double *latitude_ico = malloc(12*sizeof(double));
    double *longitude_ico = malloc(12*sizeof(double));
    int edge_vertices[NUMBER_OF_EDGES][2];
    int face_vertices[20][3];
    int face_edges[20][3];
    int face_edges_reverse[20][3];
    printf("Building icosahedron ... ");
	build_icosahedron(latitude_ico, longitude_ico, edge_vertices, face_vertices, face_edges, face_edges_reverse);
    printf("finished.\n");
    printf("Allocating memory ... ");
    double *x_unity = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *y_unity = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *z_unity = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *latitude_scalar = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *z_scalar = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *gravity_potential = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *z_vector = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *normal_distance = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *latitude_vector = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *longitude_vector = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *direction = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *gravity = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *volume = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *area = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *trsk_modified_weights = calloc(10*NUMBER_OF_VECTORS_H, sizeof(double));
    double *recov_hor_par_curl_weight = malloc(2*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_hor_ver_curl_weight = malloc(2*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_ver_0_pri_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *recov_ver_1_pri_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *recov_ver_0_curl_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *recov_ver_1_curl_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *latitude_scalar_dual = malloc(NUMBER_OF_DUAL_SCALARS_H*sizeof(double));
    double *longitude_scalar_dual = malloc(NUMBER_OF_DUAL_SCALARS_H*sizeof(double));
    double *z_scalar_dual = malloc(NUMBER_OF_DUAL_SCALARS*sizeof(double));
    double *latitude_vector_dual = malloc(NUMBER_OF_DUAL_VECTORS_PER_LAYER*sizeof(double));
    double *z_vector_dual = malloc(NUMBER_OF_DUAL_VECTORS*sizeof(double));
    double *normal_distance_dual = malloc(NUMBER_OF_DUAL_VECTORS*sizeof(double));
    double *direction_dual = malloc(NUMBER_OF_DUAL_VECTORS_H*sizeof(double));
    double *area_dual_pre = malloc(NUMBER_OF_DUAL_VECTORS*sizeof(double));
    double *area_dual = malloc((NUMBER_OF_DUAL_H_VECTORS + NUMBER_OF_H_VECTORS)*sizeof(double));
    double *f_vec = malloc((NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H)*sizeof(double));
    double *triangle_face_unit_sphere = malloc(NUMBER_OF_DUAL_VECTORS_V*sizeof(double));
    double *pent_hex_face_unity_sphere = malloc(NUMBER_OF_VECTORS_V*sizeof(double));
    double *rel_on_line_dual = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
	double *z_surface = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
	double *vertical_contravar_unit = malloc(3*NUMBER_OF_VECTORS_V*(NUMBER_OF_ORO_LAYERS + 1)*sizeof(double));
	double *e_kin_weights = malloc(14*NUMBER_OF_SCALARS*sizeof(double));
	int *e_kin_indices = malloc(14*NUMBER_OF_SCALARS*sizeof(double));
    int *to_index = malloc(NUMBER_OF_VECTORS_H*sizeof(int));
    int *from_index = malloc(NUMBER_OF_VECTORS_H*sizeof(int));
    int *trsk_modified_velocity_indices = calloc(10*NUMBER_OF_VECTORS_H, sizeof(int));
    int *trsk_modified_curl_indices = calloc(10*NUMBER_OF_VECTORS_H, sizeof(int));
    int *recov_hor_ver_pri_index = malloc(4*NUMBER_OF_VECTORS_H*sizeof(int));
    int *recov_hor_par_curl_index = malloc(2*NUMBER_OF_VECTORS_H*sizeof(int));
    int *recov_hor_ver_curl_index = malloc(2*NUMBER_OF_VECTORS_H*sizeof(int));
    int *recov_ver_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(int));
    int *adjacent_vector_indices_h = malloc(6*NUMBER_OF_SCALARS_H*sizeof(int));
    int *adjacent_vector_indices_dual_h = malloc(3*NUMBER_OF_DUAL_SCALARS_H*sizeof(int));
    int *vorticity_indices_pre = malloc(3*NUMBER_OF_DUAL_VECTORS_V*sizeof(int));
    int *vorticity_indices = malloc(4*NUMBER_OF_VECTORS_H*sizeof(int));
    int *h_curl_indices = malloc(4*NUMBER_OF_DUAL_VECTORS_H*sizeof(int));
    int *to_index_dual = malloc(NUMBER_OF_DUAL_VECTORS_H*sizeof(int));
    int *from_index_dual = malloc(NUMBER_OF_DUAL_VECTORS_H*sizeof(int));
    int *adjacent_scalar_indices_dual_h = malloc(3*NUMBER_OF_DUAL_SCALARS_H*sizeof(int));
    int *adjacent_signs_h = malloc(6*NUMBER_OF_SCALARS_H*sizeof(int));
    int *vorticity_signs_pre = malloc(3*NUMBER_OF_DUAL_VECTORS_V*sizeof(int));
    int *vorticity_signs = malloc(4*NUMBER_OF_VECTORS_H*sizeof(int));
    int *h_curl_signs = malloc(4*NUMBER_OF_DUAL_VECTORS_H*sizeof(int));
    printf("finished.\n");
    double lat_res, lon_res, base_area, base_distance, radius_0, radius_1, direction_change, x_point_0, y_point_0, z_point_0, x_point_1, y_point_1, z_point_1, x_res, y_res, z_res;
    int face_index, face_index_0, face_index_1, h_index, on_face_index, layer_index, upper_index, lower_index, coord_0_points_amount, coord_0, coord_1, triangle_on_face_index, on_edge_index, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_index, counter, primal_vector_index, dual_vector_index, number_of_triangles_per_face, edgepoint_0, edgepoint_1, edgepoint_2, points_per_edge, dual_scalar_on_face_index, base_index_old, base_index_down_triangles, base_index_up_triangles, old_triangle_on_line_index, first_face_found, edge_rel_to_face_0, edge_rel_to_face_1, edge_index, sign, small_triangle_edge_index, dump, points_upwards, points_downwards, last_triangle_bool, ncid, retval, z_surface_id, test_index;
    printf("Reading orography data ... ");
    if (ORO_ID != 0)
    {
		int ORO_FILE_LENGTH = 100;
		char *ORO_FILE_PRE = malloc((ORO_FILE_LENGTH + 1)*sizeof(char));
		sprintf(ORO_FILE_PRE, "../orography_generator/nc_files/B%d_O%d.nc", RES_ID, ORO_ID);
		ORO_FILE_LENGTH = strlen(ORO_FILE_PRE);
		free(ORO_FILE_PRE);
		char *ORO_FILE = malloc((ORO_FILE_LENGTH + 1)*sizeof(char));
		sprintf(ORO_FILE, "../orography_generator/nc_files/B%d_O%d.nc", RES_ID, ORO_ID);	    
		if ((retval = nc_open(ORO_FILE, NC_NOWRITE, &ncid)))
		    ERR(retval);
		if ((retval = nc_inq_varid(ncid, "z_surface", &z_surface_id)))
		    ERR(retval);		
		if ((retval = nc_get_var_double(ncid, z_surface_id, &z_surface[0])))
		    ERR(retval);
		if ((retval = nc_close(ncid)))
		    ERR(retval);
    }
    else
    {
    	for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
    	{
    		z_surface[i] = 0;
    	}
    }
    printf("finished.\n");
    printf("Establishing horizontal grid structure ... ");
    if (NUMBER_OF_VECTORS_H != NUMBER_OF_DUAL_VECTORS_H)
        printf("It is NUMBER_OF_VECTORS_H != NUMBER_OF_DUAL_VECTORS_H.\n");
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
	double z_oro_off = TOA*(NUMBER_OF_ORO_LAYERS + 0.0)/NUMBER_OF_LAYERS;
	double z_vertical_vector_pre[NUMBER_OF_LAYERS + 1];
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
		h_index = i - layer_index*NUMBER_OF_SCALARS_H;
		for (int j = 0; j < NUMBER_OF_LAYERS + 1; ++j)
		{
			if (j >= NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)
				z_vertical_vector_pre[j] = z_surface[h_index] + (z_oro_off - z_surface[h_index])/NUMBER_OF_ORO_LAYERS*(NUMBER_OF_LAYERS - j);
			else
				z_vertical_vector_pre[j] = TOA - (TOA - z_oro_off)/(NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)*j;
		}
		z_scalar[i] = 0.5*(z_vertical_vector_pre[layer_index] + z_vertical_vector_pre[layer_index + 1]);
        if (z_scalar[i] <= 0)
		{
            printf("z_scalar contains a non-positive value.\n");
			exit(1);
		}
    }
    double triangle_face;
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
            face_index = (i - NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
            on_face_index = i - (NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
            triangle_on_face_index = on_face_index/3;
            retval += find_coords_from_triangle_on_face_index(triangle_on_face_index, RES_ID, &coord_0, &coord_1, &coord_0_points_amount);
            dual_scalar_index = dual_scalar_on_face_index + face_index*NUMBER_OF_TRIANGLES/NUMBER_OF_BASIC_TRIANGLES;
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
            retval += find_voronoi_center_sphere(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_2], longitude_scalar[point_2], &lat_res, &lon_res);
            latitude_scalar_dual[dual_scalar_index] = lat_res;
            longitude_scalar_dual[dual_scalar_index] = lon_res;
            retval += calc_triangle_face(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_2], longitude_scalar[point_2], &triangle_face);
            triangle_face_unit_sphere[dual_scalar_index] = triangle_face;
	        retval += find_voronoi_center_sphere(latitude_scalar[point_3], longitude_scalar[point_3], latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_2], longitude_scalar[point_2], &lat_res, &lon_res);
            latitude_scalar_dual[dual_scalar_index - 1] = lat_res;
            longitude_scalar_dual[dual_scalar_index - 1] = lon_res;
            retval += calc_triangle_face(latitude_scalar[point_3], longitude_scalar[point_3], latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_2], longitude_scalar[point_2], &triangle_face);
            triangle_face_unit_sphere[dual_scalar_index - 1] = triangle_face;
            if (coord_0 == coord_0_points_amount - 1)
            {
                retval += find_voronoi_center_sphere(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_4], longitude_scalar[point_4], latitude_scalar[point_1], longitude_scalar[point_1], &lat_res, &lon_res);
                latitude_scalar_dual[dual_scalar_index + 1] = lat_res;
                longitude_scalar_dual[dual_scalar_index + 1] = lon_res;
                retval += calc_triangle_face(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_4], longitude_scalar[point_4], latitude_scalar[point_1], longitude_scalar[point_1], &triangle_face);
                triangle_face_unit_sphere[dual_scalar_index + 1] = triangle_face;
                if (coord_1 == POINTS_PER_EDGE - 1)
                {
                    retval += find_voronoi_center_sphere(latitude_scalar[point_2], longitude_scalar[point_2], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_5], longitude_scalar[point_5], &lat_res, &lon_res);
                    latitude_scalar_dual[dual_scalar_index + 2] = lat_res;
                    longitude_scalar_dual[dual_scalar_index + 2] = lon_res;
                    retval += calc_triangle_face(latitude_scalar[point_2], longitude_scalar[point_2], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_5], longitude_scalar[point_5], &triangle_face);
                    triangle_face_unit_sphere[dual_scalar_index + 2] = triangle_face;
                }
            }
        }
        retval = find_global_normal(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], &x_point_0, &y_point_0, &z_point_0);
        retval = find_global_normal(latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], &x_point_1, &y_point_1, &z_point_1);
        retval = find_between_point(x_point_0, y_point_0, z_point_0, x_point_1, y_point_1, z_point_1, 0.5, &x_res, &y_res, &z_res);
        retval = find_geos(x_res, y_res, z_res, &lat_res, &lon_res);
        latitude_vector[i] = lat_res;
        longitude_vector[i] = lon_res;
        direction[i] = find_geodetic_direction(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], 0.5);
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
    if (fabs(triangle_sum_unit_sphere/(4*M_PI) - 1) > 1e-11)
	{
        printf("Sum of faces of triangles on unit sphere does not match face of unit sphere.\n");
		exit(1);
	}
    free(x_unity);
    free(y_unity);
    free(z_unity);
    edge_rel_to_face_0 = 0;
    edge_rel_to_face_1 = 0;
    face_index_0 = 0;
    face_index_1 = 0;
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_PER_LAYER; ++i)
    {
        if (i >= NUMBER_OF_DUAL_VECTORS_H)
            latitude_vector_dual[i] = latitude_scalar_dual[i - NUMBER_OF_DUAL_VECTORS_H];
        else
        {
            if (i < NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))
            {
                edge_index = i/(POINTS_PER_EDGE + 1);
                on_edge_index = i - edge_index*(POINTS_PER_EDGE + 1);
                first_face_found = 0;
                for (int j = 0; j < NUMBER_OF_BASIC_TRIANGLES; ++j)
                {
                    if (face_edges[j][0] == edge_index || face_edges[j][1] == edge_index || face_edges[j][2] == edge_index)
                    {
                        if (first_face_found == 0)
                        {
                            face_index_0 = j;
                            first_face_found = 1;
                        }
                        else
                            face_index_1 = j;
                    }
                }
                if (face_edges[face_index_0][0] == edge_index)
                    edge_rel_to_face_0 = 0;
                if (face_edges[face_index_0][1] == edge_index)
                    edge_rel_to_face_0 = 1;
                if (face_edges[face_index_0][2] == edge_index)
                    edge_rel_to_face_0 = 2;
                if (face_edges[face_index_1][0] == edge_index)
                    edge_rel_to_face_1 = 0;
                if (face_edges[face_index_1][1] == edge_index)
                    edge_rel_to_face_1 = 1;
                if (face_edges[face_index_1][2] == edge_index)
                    edge_rel_to_face_1 = 2;
                if (edge_rel_to_face_0 == 0)
                {
                    if (face_edges_reverse[face_index_0][edge_rel_to_face_0] == 0)
                        triangle_on_face_index = 2*on_edge_index;
                    else
                        triangle_on_face_index = 2*POINTS_PER_EDGE - 2*on_edge_index;
                }
                if (edge_rel_to_face_0 == 1)
                {
                    if (face_edges_reverse[face_index_0][edge_rel_to_face_0] == 0)
                        triangle_on_face_index = -1 + (on_edge_index + 1)*(2*POINTS_PER_EDGE - on_edge_index + 1);
                    else
                        triangle_on_face_index = TRIANGLES_PER_FACE - on_edge_index*on_edge_index - 1;
                }
                if (edge_rel_to_face_0 == 2)
                {
                    if (face_edges_reverse[face_index_0][edge_rel_to_face_0] == 0)
                        triangle_on_face_index = TRIANGLES_PER_FACE - 1 - on_edge_index*(on_edge_index + 2);
                    else
                        triangle_on_face_index = on_edge_index*(2*POINTS_PER_EDGE + 2 - on_edge_index);
                }
                to_index_dual[i] = face_index_0*TRIANGLES_PER_FACE + triangle_on_face_index;
                if (edge_rel_to_face_1 == 0)
                {
                    if (face_edges_reverse[face_index_1][edge_rel_to_face_1] == 0)
                        triangle_on_face_index = 2*on_edge_index;
                    else
                        triangle_on_face_index = 2*POINTS_PER_EDGE - 2*on_edge_index;
                }
                if (edge_rel_to_face_1 == 1)
                {
                    if (face_edges_reverse[face_index_1][edge_rel_to_face_1] == 0)
                        triangle_on_face_index = -1 + (on_edge_index + 1)*(2*POINTS_PER_EDGE - on_edge_index + 1);
                    else
                        triangle_on_face_index = TRIANGLES_PER_FACE - on_edge_index*on_edge_index - 1;
                }
                if (edge_rel_to_face_1 == 2)
                {
                    if (face_edges_reverse[face_index_1][edge_rel_to_face_1] == 0)
                        triangle_on_face_index = TRIANGLES_PER_FACE - 1 - on_edge_index*(on_edge_index + 2);
                    else
                        triangle_on_face_index = on_edge_index*(2*POINTS_PER_EDGE + 2 - on_edge_index);
                }
                from_index_dual[i] = face_index_1*TRIANGLES_PER_FACE + triangle_on_face_index;
            }
            else
            {
                face_index = (i - NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
                on_face_index = i - (NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
                triangle_on_face_index = on_face_index/3;
                small_triangle_edge_index = on_face_index - 3*triangle_on_face_index;
                retval = find_coords_from_triangle_on_face_index(triangle_on_face_index, RES_ID, &coord_0, &coord_1, &coord_0_points_amount);
                if (small_triangle_edge_index == 0)
                {
                    from_index_dual[i] = face_index*TRIANGLES_PER_FACE + 2*triangle_on_face_index + coord_1;
                    to_index_dual[i] = from_index_dual[i] + 1;
                }
                if (small_triangle_edge_index == 1)
                {
                    from_index_dual[i] = face_index*TRIANGLES_PER_FACE + 2*triangle_on_face_index + 1 + coord_1;
                    to_index_dual[i] = from_index_dual[i] + 1;
                }
                if (small_triangle_edge_index == 2)
                {
                    from_index_dual[i] = face_index*TRIANGLES_PER_FACE + 2*triangle_on_face_index + 1 + coord_1;
                    to_index_dual[i] = from_index_dual[i] + 2*coord_0_points_amount;
                }
            }
            latitude_vector_dual[i] = latitude_vector[i];
            retval = find_min_dist_rel_on_line(latitude_scalar_dual[from_index_dual[i]], longitude_scalar_dual[from_index_dual[i]], latitude_scalar_dual[to_index_dual[i]], longitude_scalar_dual[to_index_dual[i]], latitude_vector_dual[i], longitude_vector[i], &rel_on_line_dual[i]);
            if (fabs(rel_on_line_dual[i] - 0.5) > 0.14)
                printf("Bisection warning.\n");
            direction_dual[i] = find_geodetic_direction(latitude_scalar_dual[from_index_dual[i]], longitude_scalar_dual[from_index_dual[i]], latitude_scalar_dual[to_index_dual[i]], longitude_scalar_dual[to_index_dual[i]], rel_on_line_dual[i]);
        }
    }
    printf("finished.\n");
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H; ++i)
    {
        if (i >= NUMBER_OF_DUAL_VECTORS_H)
        	f_vec[i] = 2*OMEGA*sin(latitude_vector[i - NUMBER_OF_DUAL_VECTORS_H]);
        else
        	f_vec[i] = 2*OMEGA*cos(latitude_vector_dual[i])*sin(direction_dual[i]);
    }
    int trouble_detected = 0;
    for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
    {
        counter = 0;
        for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
        {
            if (from_index[j] == i || to_index[j] == i)
            {
                if (from_index[j] == to_index[j])
				{
                    printf("It is from_index == to_index at some point.\n");
					exit(1);
				}
                adjacent_vector_indices_h[6*i + counter] = j;
                sign = 1;
                if (from_index[j] == i)
                    adjacent_signs_h[6*i + counter] = 1;
                if (to_index[j] == i)
                    adjacent_signs_h[6*i + counter] = -1;
                ++counter;
            }
        }
        if (counter != 6)
        {
            trouble_detected = 1;
            if (counter == 5 && i < NUMBER_OF_PENTAGONS)
                trouble_detected = 0;
        }
        if (trouble_detected == 1)
		{
            printf("Trouble detected, place 1.\n");
			exit(1);
		}
        if (i < NUMBER_OF_PENTAGONS)
        {
            adjacent_vector_indices_h[6*i + 5] = -1;
            adjacent_signs_h[6*i + 5] = 0;
        }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_SCALARS_H; ++i)
    {
    	counter = 0;
    	for (int j = 0; j < NUMBER_OF_DUAL_VECTORS_H; ++j)
    	{
    		if (from_index_dual[j] == i || to_index_dual[j] == i)
    		{
                if (from_index_dual[j] == to_index_dual[j])
				{
                    printf("It is from_index_dual == to_index_dual at some point.\n");
					exit(1);
				}
                adjacent_vector_indices_dual_h[3*i + counter] = j;
                counter++;
            }
    	}
    	if (counter != 3)
    	{
    		printf("Error in adjacent_vector_indices_dual_h creation.\n");
    		exit(1);
    	}
    }
    int number_of_edges, double_check, sign_sum_check;
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        counter = 0;
        sign_sum_check = 0;
        for (int j = 0; j < NUMBER_OF_SCALARS_H; ++j)
        {
            number_of_edges = 6;
            if (j < NUMBER_OF_PENTAGONS)
                number_of_edges = 5;
            double_check = 0;
            for (int k = 0; k < number_of_edges; ++k)
            {
                if (adjacent_vector_indices_h[6*j + k] == i)
                {
                    ++counter;
                    ++double_check;
                    sign_sum_check += adjacent_signs_h[6*j + k];
                }
            }
            if (double_check > 1)
			{
                printf("Same vector twice in adjacent_vector_indices_h of same grid cell.\n");
				exit(1);
			}
        }
        if (sign_sum_check != 0)
            printf("Problem with adjacent_signs_h.\n");
        if (counter != 2)
            printf("Problem with adjacent_vector_indices_h.\n");
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_V; ++i)
    {
        counter = 0;
        for (int j = 0; j < NUMBER_OF_DUAL_VECTORS_H; ++j)
        {
            if (from_index_dual[j] == i || to_index_dual[j] == i)
            {
                vorticity_indices_pre[3*i + counter] = j;
                sign = 1;
                if (from_index_dual[j] == i)
                {
                    find_angle_change(direction_dual[j], direction[j], &direction_change);
                    if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
                        sign = -1;
                    adjacent_scalar_indices_dual_h[3*i + counter] = to_index_dual[j];
                }
                if (to_index_dual[j] == i)
                {
                    find_angle_change(direction_dual[j], direction[j], &direction_change);
                    if (rad2deg(direction_change) > ORTH_CRITERION_DEG)
                        sign = -1;
                    adjacent_scalar_indices_dual_h[3*i + counter] = from_index_dual[j];
                }
                vorticity_signs_pre[3*i + counter] = sign;
                ++counter;
            }
        }
        if (counter != 3)
		{
            printf("Trouble detected, place 0.\n");
			exit(1);
		}
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        find_angle_change(direction[i], direction_dual[i], &direction_change);
        if (fabs(rad2deg(direction_change)) < ORTH_CRITERION_DEG || fabs(rad2deg(direction_change)) > 90 + (90 - ORTH_CRITERION_DEG))
		{
            printf("Grid non-orthogonal.\n");
			exit(1);
		}
    }
    int check_0, check_1, check_2;
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
	int dual_scalar_h_index;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_VECTORS_V)
        {
            gravity[i] = 0;
            z_vector[i] = 0.5*(z_scalar[layer_index*NUMBER_OF_SCALARS_H + to_index[h_index - NUMBER_OF_VECTORS_V]] + z_scalar[layer_index*NUMBER_OF_SCALARS_H + from_index[h_index - NUMBER_OF_VECTORS_V]]);
            if (z_vector[i] <= 0)
			{
                printf("z_vector contains a non-positive value at a horizontal grid point.\n");
				exit(1);
			}
            normal_distance[i] = calculate_distance_h(latitude_scalar[from_index[h_index - NUMBER_OF_VECTORS_V]], longitude_scalar[from_index[h_index - NUMBER_OF_VECTORS_V]], latitude_scalar[to_index[h_index - NUMBER_OF_VECTORS_V]], longitude_scalar[to_index[h_index - NUMBER_OF_VECTORS_V]], RADIUS + z_vector[i]);
        }
        else
        {
            upper_index = h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H;
            lower_index = h_index + layer_index*NUMBER_OF_SCALARS_H;
            if (layer_index == 0)
			{
                normal_distance[i] = TOA - z_scalar[lower_index];
            	z_vector[i] = TOA;
			}
            else if (layer_index == NUMBER_OF_LAYERS)
			{
                normal_distance[i] = z_scalar[upper_index] - z_surface[h_index];
				z_vector[i] = z_surface[h_index];
			}
            else
			{
                normal_distance[i] = z_scalar[upper_index] - z_scalar[lower_index];
				z_vector[i] = z_scalar[lower_index] + 0.5*normal_distance[i];
			}
            if (z_vector[i] < z_surface[h_index])
			{
                printf("z_vector lays below surface.\n");
				exit(1);
			}
            area[i] = pent_hex_face_unity_sphere[h_index]*pow(RADIUS + z_vector[i], 2);
            gravity[i] = -GRAVITY_MEAN_SFC_ABS*pow(RADIUS, 2)/pow(RADIUS + z_vector[i], 2);
        }
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        if (normal_distance[i] <= 0)
		{
            printf("normal_distance contains a non-positive value.\n");
			exit(1);
		}
    }
	int index_vector_for_dual_scalar_z[3];
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_VECTORS_V)
        {
            if (h_index - NUMBER_OF_VECTORS_V >= NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))
            {
                face_index = (h_index - NUMBER_OF_VECTORS_V - NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
                on_face_index = h_index - NUMBER_OF_VECTORS_V - (NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
                triangle_on_face_index = on_face_index/3;
                retval = find_coords_from_triangle_on_face_index(triangle_on_face_index, RES_ID, &coord_0, &coord_1, &coord_0_points_amount);
				dual_scalar_h_index = face_index*TRIANGLES_PER_FACE + 1 + 2*triangle_on_face_index + coord_1;
                dual_scalar_index = layer_index*NUMBER_OF_DUAL_SCALARS_H + dual_scalar_h_index;
				retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index, index_vector_for_dual_scalar_z);
                z_scalar_dual[dual_scalar_index] = 1.0/3*(z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
				retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index - 1, index_vector_for_dual_scalar_z);
                z_scalar_dual[dual_scalar_index - 1] = 1.0/3*(z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                if (layer_index == NUMBER_OF_LAYERS - 1)
                {
					retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index, index_vector_for_dual_scalar_z);
                	z_scalar_dual[dual_scalar_index + NUMBER_OF_DUAL_SCALARS_H] = 1.0/3*(z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
					retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index - 1, index_vector_for_dual_scalar_z);
                	z_scalar_dual[dual_scalar_index + NUMBER_OF_DUAL_SCALARS_H - 1] = 1.0/3*(z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                }
                if (coord_0 == coord_0_points_amount - 1)
                {
					retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index + 1, index_vector_for_dual_scalar_z);
                	z_scalar_dual[dual_scalar_index + 1] = 1.0/3*(z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                    if (layer_index == NUMBER_OF_LAYERS - 1)
            			z_scalar_dual[dual_scalar_index + 1 + NUMBER_OF_DUAL_SCALARS_H] = 1.0/3*(z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                    if (coord_1 == POINTS_PER_EDGE - 1)
                    {
						retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, dual_scalar_h_index + 2, index_vector_for_dual_scalar_z);
                		z_scalar_dual[dual_scalar_index + 2] = 1.0/3*(z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[layer_index*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                        if (layer_index == NUMBER_OF_LAYERS - 1)
            				z_scalar_dual[dual_scalar_index + 2 + NUMBER_OF_DUAL_SCALARS_H] = 1.0/3*(z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[0]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[1]] + z_vector[(layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER + index_vector_for_dual_scalar_z[2]]);
                    }
                }
            }
        }
    }
	int min_oro_index = find_min_index(z_surface, NUMBER_OF_SCALARS_H);
	double min_oro = z_surface[min_oro_index];
	for (int i = 0; i < NUMBER_OF_DUAL_SCALARS; ++i)
	{
		if (z_scalar_dual[i] > TOA)
		{
			printf("z_scalar_dual has a point above top of atmosphere.\n");
			exit(1);
		}
		if (z_scalar_dual[i] < min_oro)
		{
			printf("z_scalar_dual has a point below minimum of orography.\n");
			exit(1);
		}
	}
    free(pent_hex_face_unity_sphere);
	double check_sum;
	for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
	{
		check_sum = 0;
		for (int j = 0; j < NUMBER_OF_LEVELS; ++j)
		{
			check_sum += normal_distance[i + j*NUMBER_OF_VECTORS_PER_LAYER];
		}
		if (fabs(check_sum/(TOA - z_surface[i]) - 1) > 1e-15)
		{
			printf("Problem 0 with vertical grid structure.\n");
			exit(1);
		}
	}
    double volume_sum, volume_sum_ideal;
    volume_sum = 0;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
        h_index = i - layer_index*NUMBER_OF_SCALARS_H;
        base_area = area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        radius_0 = RADIUS + z_vector[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        radius_1 = RADIUS + z_vector[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER];
        volume[i] = find_volume(base_area, radius_0, radius_1);
        volume_sum += volume[i];
    }
    volume_sum_ideal = 0;
    for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
    	volume_sum_ideal += find_volume(area[NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V + i], RADIUS + z_vector[NUMBER_OF_VECTORS- NUMBER_OF_VECTORS_V + i], RADIUS + TOA);
    if (fabs(volume_sum/volume_sum_ideal - 1) > 1e-12)
	{
        printf("Sum of volumes of grid boxes does not match volume of entire atmosphere.\n");
		exit(1);
	}
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        if (volume[i] <= 0)
		{
            printf("volume contains a non-positive value.\n");
			exit(1);
		}
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_DUAL_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_DUAL_VECTORS_H)
        {
            upper_index = h_index - NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_SCALARS_H;
            lower_index = h_index - NUMBER_OF_DUAL_VECTORS_H + (layer_index + 1)*NUMBER_OF_DUAL_SCALARS_H;
            normal_distance_dual[i] = z_scalar_dual[upper_index] - z_scalar_dual[lower_index];
			z_vector_dual[i] = z_scalar_dual[lower_index] + 0.5*normal_distance_dual[i];
            area_dual_pre[i] = pow(RADIUS + z_vector_dual[i], 2)*triangle_face_unit_sphere[h_index - NUMBER_OF_DUAL_VECTORS_H];
        }
        else
        {
			if (layer_index == 0)
				z_vector_dual[i] = TOA;
			else if (layer_index == NUMBER_OF_LAYERS)
				z_vector_dual[i] = 0.5*(z_surface[from_index[h_index]] + z_surface[to_index[h_index]]);
			else
				z_vector_dual[i] = 0.5*(z_vector[NUMBER_OF_VECTORS_V + h_index + (layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER] + z_vector[NUMBER_OF_VECTORS_V + h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER]);
            if (layer_index == 0)
                radius_1 = RADIUS + z_vector_dual[i];
            else
                radius_1 = RADIUS + z_vector[NUMBER_OF_VECTORS_V + h_index + (layer_index - 1)*NUMBER_OF_VECTORS_PER_LAYER];
            if (layer_index == NUMBER_OF_LAYERS)
                radius_0 = RADIUS + z_vector_dual[i];
            else
                radius_0 = RADIUS + z_vector[NUMBER_OF_VECTORS_V + h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER];
            primal_vector_index = (NUMBER_OF_LAYERS - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + h_index;
            base_distance = normal_distance[primal_vector_index]/(RADIUS + z_vector[primal_vector_index])*radius_0;
            area_dual_pre[i] = calculate_vertical_face(base_distance, radius_0, radius_1);
            normal_distance_dual[i] = calculate_distance_h(latitude_scalar_dual[from_index_dual[h_index]], longitude_scalar_dual[from_index_dual[h_index]], latitude_scalar_dual[to_index_dual[h_index]], longitude_scalar_dual[to_index_dual[h_index]], RADIUS + z_vector_dual[i]);
        }
    }
    free(triangle_face_unit_sphere);
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; ++i)
    {
        if (area_dual_pre[i] <= 0)
		{
            printf("area_dual_pre contains a non-positive value.\n");
			exit(1);
		}
        if (normal_distance_dual[i] <= 0)
		{
            printf("normal_distance_dual contains a non-positive value.\n");
			exit(1);
		}
    }
	for (int i = 0; i < NUMBER_OF_DUAL_SCALARS_H; ++i)
	{
		check_sum = 0;
		for (int j = 0; j < NUMBER_OF_LAYERS; ++j)
		{
			check_sum += normal_distance_dual[i + NUMBER_OF_DUAL_VECTORS_H + j*NUMBER_OF_DUAL_VECTORS_PER_LAYER];
		}
		retval = find_v_vector_indices_for_dual_scalar_z(from_index, to_index, vorticity_indices_pre, i, index_vector_for_dual_scalar_z);
		if (fabs(check_sum/(TOA - 1.0/3*(z_surface[index_vector_for_dual_scalar_z[0]] + z_surface[index_vector_for_dual_scalar_z[1]] + z_surface[index_vector_for_dual_scalar_z[2]])) - 1) > 1e-15)
		{
			printf("Problem 1 with vertical grid structure.\n");
			exit(1);
		}
	}
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_VECTORS_V)
        {
            dual_vector_index = NUMBER_OF_LAYERS*NUMBER_OF_DUAL_VECTORS_PER_LAYER + h_index - NUMBER_OF_VECTORS_V;
            radius_1 = RADIUS + z_vector_dual[h_index - NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER];
            radius_0 = RADIUS + z_vector_dual[h_index - NUMBER_OF_VECTORS_V + (layer_index + 1)*NUMBER_OF_DUAL_VECTORS_PER_LAYER];
            base_distance = normal_distance_dual[dual_vector_index]/(RADIUS + z_vector_dual[dual_vector_index])*radius_0;
            area[i] = calculate_vertical_face(base_distance, radius_0, radius_1);
        }
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        if (area[i] <= 0)
		{
            printf("It is area <= 0 at some point, position 0.\n");
			exit(1);
		}
    }
    double cell_base_area;
    double alpha_1 = 3.0/4;
    alpha_1 = 1;
    double alpha_2 = 1 - alpha_1;
    double triangle_face_0, triangle_face_1;
    double z_value_0, z_value_1, z_value_2, z_triangle_mean, partial_volume, kite_area_0, kite_area_1, kite_area, triangle_area_0, triangle_area_1, triangle_area_0_for_kite, triangle_area_1_for_kite, triangle_area_for_kite;
    int second_index, index_found, neighbouring_cell_h_index, vertex_index;
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
				retval = calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + j]], longitude_vector[adjacent_vector_indices_h[6*h_index + j]], latitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_face_0);
				if (retval != 0)
				{
					printf("Error in e_kin calculation, position 0.\n");
					exit(1);
				}
				z_value_0 = z_scalar[i];
				z_value_1 = z_vector[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				z_value_2 = z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + from_index_dual[adjacent_vector_indices_h[6*h_index + j]]];
				z_triangle_mean = 1.0/3*(z_value_0 + z_value_1 + z_value_2);
				triangle_face_0 = triangle_face_0*pow(RADIUS + z_triangle_mean, 2);
				retval = calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + j]], longitude_vector[adjacent_vector_indices_h[6*h_index + j]], latitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_face_1);
				if (retval != 0)
				{
					printf("Error in e_kin calculation, position 1.\n");
					exit(1);
				}
				z_value_0 = z_scalar[i];
				z_value_1 = z_vector[NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				z_value_2 = z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + to_index_dual[adjacent_vector_indices_h[6*h_index + j]]];
				z_triangle_mean = 1.0/3*(z_value_0 + z_value_1 + z_value_2);
				triangle_face_1 = triangle_face_1*pow(RADIUS + z_triangle_mean, 2);
				retval = calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + j]], longitude_vector[adjacent_vector_indices_h[6*h_index + j]], latitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area_0);
				if (retval != 0)
				{
					printf("Error in e_kin calculation, position 2.\n");
					exit(1);
				}
				second_index = -1;
				for (int k = 0; k < number_of_edges; ++k)
				{
					if (k != j && (to_index_dual[adjacent_vector_indices_h[6*h_index + k]] == from_index_dual[adjacent_vector_indices_h[6*h_index + j]] || from_index_dual[adjacent_vector_indices_h[6*h_index + k]] == from_index_dual[adjacent_vector_indices_h[6*h_index + j]]))
						second_index = k;
				}
				if (second_index == -1)
				{
					printf("Error with second_index in kinetic energy calculation, position 0.\n");
					exit(1);
				}
				retval = calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + second_index]], longitude_vector[adjacent_vector_indices_h[6*h_index + second_index]], latitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area_1);
				if (retval != 0)
				{
					printf("Error in e_kin calculation, position 3.\n");
					exit(1);
				}
				kite_area_0 = triangle_area_0 + triangle_area_1;
				retval = calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + j]], longitude_vector[adjacent_vector_indices_h[6*h_index + j]], latitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area_0);
				if (retval != 0)
				{
					printf("Error in e_kin calculation, position 4.\n");
					exit(1);
				}
				second_index = -1;
				for (int k = 0; k < number_of_edges; ++k)
				{
					if (k != j && (to_index_dual[adjacent_vector_indices_h[6*h_index + k]] == to_index_dual[adjacent_vector_indices_h[6*h_index + j]] || from_index_dual[adjacent_vector_indices_h[6*h_index + k]] == to_index_dual[adjacent_vector_indices_h[6*h_index + j]]))
						second_index = k;
				}
				if (second_index == -1)
				{
					printf("Error with second_index in kinetic energy calculation, position 1.\n");
					exit(1);
				}
				retval = calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + second_index]], longitude_vector[adjacent_vector_indices_h[6*h_index + second_index]], latitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area_1);
				if (retval != 0)
				{
					printf("Error in e_kin calculation, position 5.\n");
					exit(1);
				}
				kite_area_1 = triangle_area_0 + triangle_area_1;
				kite_area_0 = kite_area_0*pow(RADIUS + z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], 2);
				kite_area_1 = kite_area_1*pow(RADIUS + z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], 2);
				neighbouring_cell_h_index = from_index[adjacent_vector_indices_h[6*h_index + j]];
				if (neighbouring_cell_h_index == h_index)
					neighbouring_cell_h_index = to_index[adjacent_vector_indices_h[6*h_index + j]];
				retval = calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_scalar[neighbouring_cell_h_index], longitude_scalar[neighbouring_cell_h_index], latitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area_0_for_kite);
				if (retval != 0)
				{
					printf("Error in e_kin calculation, position 8.\n");
					exit(1);
				}
				retval = calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_scalar[neighbouring_cell_h_index], longitude_scalar[neighbouring_cell_h_index], latitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], longitude_scalar_dual[to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], &triangle_area_1_for_kite);
				if (retval != 0)
				{
					printf("Error in e_kin calculation, position 9.\n");
					exit(1);
				}
				triangle_area_0_for_kite = triangle_area_0_for_kite*pow(z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + from_index_dual[adjacent_vector_indices_h[6*h_index + j]]], 2);
				triangle_area_1_for_kite = triangle_area_1_for_kite*pow(z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + to_index_dual[adjacent_vector_indices_h[6*h_index + j]]], 2);
				e_kin_weights[14*i + j] = alpha_1*(triangle_face_0 + triangle_face_1)/base_area + alpha_2*(kite_area_0/base_area*triangle_area_0_for_kite/area_dual_pre[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + from_index_dual[adjacent_vector_indices_h[6*h_index + j]]] + kite_area_1/base_area*triangle_area_1_for_kite/area_dual_pre[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + to_index_dual[adjacent_vector_indices_h[6*h_index + j]]]);
			}
			else
			{
				e_kin_indices[14*i + j] = 0;
				e_kin_weights[14*i + j] = 0;
			}
		}
    	for (int j = 6; j < 12; ++j)
    	{
    		vertex_index = from_index_dual[adjacent_vector_indices_h[6*h_index + j - 6]];
    		if (j < 11 || h_index >= NUMBER_OF_PENTAGONS)
    		{
				second_index = -1;
				for (int k = 0; k < number_of_edges; ++k)
				{
					if (k != j && (to_index_dual[adjacent_vector_indices_h[6*h_index + k]] == vertex_index || from_index_dual[adjacent_vector_indices_h[6*h_index + k]] == vertex_index))
						second_index = k;
				}
				if (second_index == -1)
				{
					printf("Error with second_index in kinetic energy calculation, position 2.\n");
					exit(1);
				}
				index_found = 0;
				if (adjacent_vector_indices_dual_h[3*vertex_index + 0] == adjacent_vector_indices_h[6*h_index + j - 6] || adjacent_vector_indices_dual_h[3*vertex_index + 0] == adjacent_vector_indices_h[6*h_index + second_index])
				{
					index_found = 1;
				}
				if (adjacent_vector_indices_dual_h[3*vertex_index + 1] == adjacent_vector_indices_h[6*h_index + j - 6] || adjacent_vector_indices_dual_h[3*vertex_index + 1] == adjacent_vector_indices_h[6*h_index + second_index])
					index_found = 2;
				e_kin_indices[14*i + j] = NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_VECTORS_PER_LAYER + adjacent_vector_indices_dual_h[3*vertex_index + index_found];
				retval = calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + j - 6]], longitude_vector[adjacent_vector_indices_h[6*h_index + j - 6]], latitude_scalar_dual[vertex_index], longitude_scalar_dual[vertex_index], &triangle_area_0);
				if (retval != 0)
				{
					printf("Error in e_kin calculation, position 6.\n");
					exit(1);
				}
				retval = calc_triangle_face(latitude_scalar[h_index], longitude_scalar[h_index], latitude_vector[adjacent_vector_indices_h[6*h_index + second_index]], longitude_vector[adjacent_vector_indices_h[6*h_index + second_index]], latitude_scalar_dual[vertex_index], longitude_scalar_dual[vertex_index], &triangle_area_1);
				if (retval != 0)
				{
					printf("Error in e_kin calculation, position 7.\n");
					exit(1);
				}
				kite_area = triangle_area_0 + triangle_area_1;
				kite_area = kite_area*pow(RADIUS + z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + vertex_index], 2);
				retval = calc_triangle_face(latitude_scalar[to_index[adjacent_vector_indices_dual_h[3*vertex_index + index_found]]], longitude_scalar[to_index[adjacent_vector_indices_dual_h[3*vertex_index + index_found]]], latitude_scalar[from_index[adjacent_vector_indices_dual_h[3*vertex_index + index_found]]], longitude_scalar[from_index[adjacent_vector_indices_dual_h[3*vertex_index + index_found]]], latitude_scalar_dual[vertex_index], longitude_scalar_dual[vertex_index], &triangle_area_for_kite);
				if (retval != 0)
				{
					printf("Error in e_kin calculation, position 10.\n");
					exit(1);
				}
				triangle_area_for_kite = triangle_area_for_kite*pow(z_vector_dual[NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + vertex_index], 2);
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
    printf("Checking kinetic energy weights sum ... ");
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
		if (fabs(e_kin_weights_sum - 1.5) > 0.000001)
		{
			printf("%lf\n", e_kin_weights_sum);
			printf("Error in e_kin_weights, position 1.\n");
			exit(1);
		}
    }
    printf("passed.\n");
    double x_value, y_value, z_value, normal_vector_global[3], local_basis_vector[3], local_x, local_y, local_z, tilting_angle, delta_x, delta_z, abs_value;
    for (int i = 0; i < NUMBER_OF_ORO_LAYERS + 1; ++i)
    {
    	for (int j = 0; j < NUMBER_OF_VECTORS_V; ++j)
    	{
    		if (j < NUMBER_OF_PENTAGONS)
    		{
    			number_of_edges = 5;
			}
    		else
    		{
    			number_of_edges = 6;
			}
			normal_vector_global[0] = 0;
			normal_vector_global[1] = 0;
			normal_vector_global[2] = 0;
    		for (int k = 0; k < number_of_edges; ++k)
    		{
    			layer_index = i + NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS;
    			dual_vector_index = layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*j + k];
    			delta_x = normal_distance_dual[dual_vector_index];
    			delta_z = z_scalar_dual[layer_index*NUMBER_OF_DUAL_SCALARS_H + to_index_dual[adjacent_vector_indices_h[6*j + k]]] - z_scalar_dual[layer_index*NUMBER_OF_DUAL_SCALARS_H + from_index_dual[adjacent_vector_indices_h[6*j + k]]];
    			tilting_angle = atan(delta_z/delta_x);
    			local_x = -sin(tilting_angle)*cos(direction_dual[adjacent_vector_indices_h[6*j + k]]);
    			local_y = -sin(tilting_angle)*sin(direction_dual[adjacent_vector_indices_h[6*j + k]]);
    			local_z = cos(tilting_angle);	
    			abs_value = sqrt(pow(local_x, 2) + pow(local_y, 2) + pow(local_z, 2));
    			local_x = local_x/abs_value;
    			local_y = local_y/abs_value;
    			local_z = local_z/abs_value;
    			calc_local_i(latitude_scalar[j], longitude_scalar[j], local_basis_vector);
    			normal_vector_global[0] += local_x*local_basis_vector[0];
    			normal_vector_global[1] += local_x*local_basis_vector[1];
    			normal_vector_global[2] += local_x*local_basis_vector[2];
    			calc_local_j(latitude_scalar[j], longitude_scalar[j], local_basis_vector);
    			normal_vector_global[0] += local_y*local_basis_vector[0];
    			normal_vector_global[1] += local_y*local_basis_vector[1];
    			normal_vector_global[2] += local_y*local_basis_vector[2];
    			calc_local_k(latitude_scalar[j], longitude_scalar[j], local_basis_vector);
    			normal_vector_global[0] += local_z*local_basis_vector[0];
    			normal_vector_global[1] += local_z*local_basis_vector[1];
    			normal_vector_global[2] += local_z*local_basis_vector[2];
    		}
    		abs_value = sqrt(scalar_product_elementary(normal_vector_global, normal_vector_global));
    		normal_vector_global[0] = normal_vector_global[0]/abs_value;
    		normal_vector_global[1] = normal_vector_global[1]/abs_value;
    		normal_vector_global[2] = normal_vector_global[2]/abs_value;
    		calc_local_i(latitude_scalar[j], longitude_scalar[j], local_basis_vector);
    		x_value = scalar_product_elementary(normal_vector_global, local_basis_vector);
    		calc_local_j(latitude_scalar[j], longitude_scalar[j], local_basis_vector);
    		y_value = scalar_product_elementary(normal_vector_global, local_basis_vector);
    		calc_local_k(latitude_scalar[j], longitude_scalar[j], local_basis_vector);
    		z_value = scalar_product_elementary(normal_vector_global, local_basis_vector);
    		for (int k = 0; k < 3; ++k)
    		{
    			if (k == 0)
    			{
    				vertical_contravar_unit[i*3*NUMBER_OF_VECTORS_V + 3*j + k] = x_value;
				}
    			if (k == 1)
    			{
    				vertical_contravar_unit[i*3*NUMBER_OF_VECTORS_V + 3*j + k] = y_value;
				}
    			if (k == 2)
    			{
    				vertical_contravar_unit[i*3*NUMBER_OF_VECTORS_V + 3*j + k] = z_value;
				}
				if (fabs(vertical_contravar_unit[i*3*NUMBER_OF_VECTORS_V + 3*j + k]) > 1.000001)
				{
				    printf("fabs(vertical_contravar_unit) > 1 at some point.\n");
					exit(1);
				}
			}
    	}
    }
    double area_rescale_factor, dz, dh;
	double area_rescale_factor_warning_begin = 0.1;
	int contravar_unit_vector_index;
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
    	layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
    	h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
    	if (h_index < NUMBER_OF_VECTORS_V && layer_index >= NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS)
    	{
    		contravar_unit_vector_index = h_index + (layer_index - (NUMBER_OF_LAYERS - NUMBER_OF_ORO_LAYERS))*NUMBER_OF_VECTORS_V;
    		dz = vertical_contravar_unit[3*contravar_unit_vector_index + 2];
    		if (dz <= 0)
    			printf("Error in horizontal face orography rescale calculation, position 0.\n");
    		dh = sqrt(pow(vertical_contravar_unit[3*contravar_unit_vector_index + 0], 2) + pow(vertical_contravar_unit[3*contravar_unit_vector_index + 1], 2));
    		tilting_angle = atan(dh/dz);
    		area_rescale_factor = sqrt(1 + pow(tan(tilting_angle), 2));
    		area[i] = area_rescale_factor*area[i];
    		if (fabs(area_rescale_factor - 1) > area_rescale_factor_warning_begin)
    			printf("It is area_rescale_factor > area_rescale_factor_warning_begin at some point.\n");
    		if (area[i] <= 0)
    		{
    			printf("It is area <= 0 at some point, position 1.\n");
    			exit(1);
    		}
    	}
    }
    /* int bool_0, bool_1, first_found;
    int cell_0_for_cross, cell_1_for_cross;
	uncomment for ald Coriolis reconstruction	
	*/    
	int *face_of_cell_indices = malloc(2*sizeof(int));
	int offset, sign_0, sign_1;
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
		double triangle_0, triangle_1, check_sum;
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
						int vertex_indices[5];
						int edge_indices[5];
						int indices_resorted[5];
						int vertex_indices_resorted[5];
						double latitude_vertices[5];
						double longitude_vertices[5];
						double latitude_edges[5];
						double longitude_edges[5];
						double vector_of_areas[5];
						for (int l = 0; l < 5; ++l)
							vertex_indices[l] = -1;
						counter = 0;
						for (int l = 0; l < 5; ++l)
						{
							vertex_index_candidate_0 = from_index_dual[adjacent_vector_indices_h[6*from_index[i] + l]];
							vertex_index_candidate_1 = to_index_dual[adjacent_vector_indices_h[6*from_index[i] + l]];
							retval = in_bool_calculator(vertex_indices, 5, vertex_index_candidate_0, &check_result);						
							if (check_result == 0)
							{
								vertex_indices[counter] = vertex_index_candidate_0;
								latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
								longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
								++counter;
							}
							retval = in_bool_calculator(vertex_indices, 5, vertex_index_candidate_1, &check_result);						
							if (check_result == 0)
							{
								vertex_indices[counter] = vertex_index_candidate_1;
								latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
								longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
								++counter;
							}
						}
						if (counter != 5)
						{
							printf("Problem 13 in TRSK implementation detected.\n");
							exit(1);
						}
						for (int l = 0; l < 5; ++l)
						{
							if (vertex_indices[l] < 0 || vertex_indices[l] >= NUMBER_OF_DUAL_SCALARS_H)
							{
								printf("Problem 3 in TRSK implementation detected.\n");
								exit(1);
							}
							for (int m = l + 1; m < 5; ++m)
							{
								if (vertex_indices[l] == vertex_indices[m])
								{
									printf("Problem 17 in TRSK implementation detected.\n");
									exit(1);
								}
							}
						}
						retval = sort_edge_indices(latitude_vertices, longitude_vertices, 5, indices_resorted);
						for (int l = 0; l < 5; ++l)
						{
							for (int m = l + 1; m < 5; ++m)
							{
								if (indices_resorted[l] == indices_resorted[m] || indices_resorted[l] < 0 || indices_resorted[l] > 4)
								{
									printf("Problem 25 in TRSK implementation detected.\n");
									exit(1);
								}
							}
						}
						for (int l = 0; l < 5; ++l)
							vertex_indices_resorted[l] = vertex_indices[indices_resorted[l]];
						for (int l = 0; l < 5; ++l)
						{
							for (int m = 0; m < 5; ++m)
							{
								if ((from_index_dual[adjacent_vector_indices_h[6*from_index[i] + m]] == vertex_indices_resorted[l] && to_index_dual[adjacent_vector_indices_h[6*from_index[i] + m]] == vertex_indices_resorted[(l + 1)%5]) || (to_index_dual[adjacent_vector_indices_h[6*from_index[i] + m]] == vertex_indices_resorted[l] && from_index_dual[adjacent_vector_indices_h[6*from_index[i] + m]] == vertex_indices_resorted[(l + 1)%5]))
									edge_indices[l] = adjacent_vector_indices_h[6*from_index[i] + m];
							}
						}
						for (int l = 0; l < 5; ++l)
						{
							if (edge_indices[l] < 0 || edge_indices[l] >= NUMBER_OF_VECTORS_H)
							{
								printf("Problem 7 in TRSK implementation detected.\n");
								exit(1);
							}
							for (int m = l + 1; m < 5; ++m)
							{
								if (edge_indices[l] == edge_indices[m])
								{
									printf("Problem 18 in TRSK implementation detected.\n");
									exit(1);
								}
							}
						}
						for (int l = 0; l < 5; ++l)
						{
							latitude_edges[l] = latitude_vector[edge_indices[l]];
							longitude_edges[l] = longitude_vector[edge_indices[l]];
						}
						check_sum = 0;
						for (int l = 0; l < 5; ++l)	
						{
							if (l == 0)
								retval = calc_triangle_face(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[4], longitude_edges[4], &triangle_0);
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
						for (int l = 0; l < 5; ++l)
						{
							if (edge_indices[l] == i)
								last_index = l;
							if (edge_indices[l] == trsk_modified_velocity_indices[10*i + k])
								first_index = (l + 1)%5;
						}
						if (k == 4)
							sum_of_weights = 0;
						else
							double_sum_gen(vector_of_areas, 5, first_index, last_index, &sum_of_weights);
						if (sum_of_weights < 0 || sum_of_weights/area[from_index[i]] > 1)
						{
							printf("Problem 34 in TRSK implementation detected\n");
							exit(1);
						}
					}
					else
					{
						int vertex_indices[6];
						int edge_indices[6];
						int indices_resorted[6];
						int vertex_indices_resorted[6];
						double latitude_vertices[6];
						double longitude_vertices[6];
						double latitude_edges[6];
						double longitude_edges[6];
						double vector_of_areas[6];
						for (int l = 0; l < 6; ++l)
							vertex_indices[l] = -1;
						counter = 0;
						for (int l = 0; l < 6; ++l)
						{
							vertex_index_candidate_0 = from_index_dual[adjacent_vector_indices_h[6*from_index[i] + l]];
							vertex_index_candidate_1 = to_index_dual[adjacent_vector_indices_h[6*from_index[i] + l]];
							retval = in_bool_calculator(vertex_indices, 6, vertex_index_candidate_0, &check_result);						
							if (check_result == 0)
							{
								vertex_indices[counter] = vertex_index_candidate_0;
								latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
								longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
								++counter;
							}
							retval = in_bool_calculator(vertex_indices, 6, vertex_index_candidate_1, &check_result);						
							if (check_result == 0)
							{
								vertex_indices[counter] = vertex_index_candidate_1;
								latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
								longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
								++counter;
							}
						}
						if (counter != 6)
						{
							printf("Problem 14 in TRSK implementation detected.\n");
							exit(1);
						}
						for (int l = 0; l < 6; ++l)
						{
							if (vertex_indices[l] < 0 || vertex_indices[l] >= NUMBER_OF_DUAL_SCALARS_H)
							{
								printf("Problem 4 in TRSK implementation detected.\n");
								exit(1);
							}
							for (int m = l + 1; m < 6; ++m)
							{
								if (vertex_indices[l] == vertex_indices[m])
								{
									printf("Problem 22 in TRSK implementation detected.\n");
									exit(1);
								}
							}
						}
						retval = sort_edge_indices(latitude_vertices, longitude_vertices, 6, indices_resorted);
						for (int l = 0; l < 6; ++l)
						{
							for (int m = l + 1; m < 6; ++m)
							{
								if (indices_resorted[l] == indices_resorted[m] || indices_resorted[l] < 0 || indices_resorted[l] > 5)
								{
									printf("Problem 26 in TRSK implementation detected.\n");
									exit(1);
								}
							}
						}
						for (int l = 0; l < 6; ++l)
							vertex_indices_resorted[l] = vertex_indices[indices_resorted[l]];
						for (int l = 0; l < 6; ++l)
						{
							for (int m = 0; m < 6; ++m)
							{
								if ((from_index_dual[adjacent_vector_indices_h[6*from_index[i] + m]] == vertex_indices_resorted[l] && to_index_dual[adjacent_vector_indices_h[6*from_index[i] + m]] == vertex_indices_resorted[(l + 1)%6]) || (to_index_dual[adjacent_vector_indices_h[6*from_index[i] + m]] == vertex_indices_resorted[l] && from_index_dual[adjacent_vector_indices_h[6*from_index[i] + m]] == vertex_indices_resorted[(l + 1)%6]))
									edge_indices[l] = adjacent_vector_indices_h[6*from_index[i] + m];
							}
						}
						for (int l = 0; l < 6; ++l)
						{
							if (edge_indices[l] < 0 || edge_indices[l] >= NUMBER_OF_VECTORS_H)
							{
								printf("Problem 8 in TRSK implementation detected.\n");
								exit(1);
							}
							for (int m = l + 1; m < 6; ++m)
							{
								if (edge_indices[l] == edge_indices[m])
								{
									printf("Problem 19 in TRSK implementation detected.\n");
									exit(1);
								}
							}
						}
						for (int l = 0; l < 6; ++l)
						{
							latitude_edges[l] = latitude_vector[edge_indices[l]];
							longitude_edges[l] = longitude_vector[edge_indices[l]];
						}
						check_sum = 0;
						for (int l = 0; l < 6; ++l)	
						{
							if (l == 0)
								retval = calc_triangle_face(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[5], longitude_edges[5], &triangle_0);
							else
								retval = calc_triangle_face(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[l - 1], longitude_edges[l - 1], &triangle_0);
							retval = calc_triangle_face(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[l], longitude_edges[l], &triangle_1);
							vector_of_areas[l] = pow(RADIUS + z_scalar[from_index[i]], 2)*(triangle_0 + triangle_1);
							check_sum += vector_of_areas[l];
						}
						if (fabs(check_sum/area[from_index[i]] - 1) > 0.001)
						{
							printf("Problem 31 in TRSK implementation detected.\n");
							exit(1);
						}
						for (int l = 0; l < 6; ++l)
						{
							if (edge_indices[l] == i)
								last_index = l;
							if (edge_indices[l] == trsk_modified_velocity_indices[10*i + k])
								first_index = (l + 1)%6;
						}
						double_sum_gen(vector_of_areas, 6, first_index, last_index, &sum_of_weights);
						if (sum_of_weights < 0 || sum_of_weights/area[from_index[i]] > 1)
						{
							printf("Problem 35 in TRSK implementation detected\n");
							exit(1);
						}
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
						int vertex_indices[5];
						int edge_indices[5];
						int indices_resorted[5];
						int vertex_indices_resorted[5];
						double latitude_vertices[5];
						double longitude_vertices[5];
						double latitude_edges[5];
						double longitude_edges[5];
						double vector_of_areas[5];
						for (int l = 0; l < 5; ++l)
							vertex_indices[l] = -1;
						counter = 0;
						for (int l = 0; l < 5; ++l)
						{
							vertex_index_candidate_0 = from_index_dual[adjacent_vector_indices_h[6*to_index[i] + l]];
							vertex_index_candidate_1 = to_index_dual[adjacent_vector_indices_h[6*to_index[i] + l]];
							retval = in_bool_calculator(vertex_indices, 5, vertex_index_candidate_0, &check_result);						
							if (check_result == 0)
							{
								vertex_indices[counter] = vertex_index_candidate_0;
								latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
								longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
								++counter;
							}
							retval = in_bool_calculator(vertex_indices, 5, vertex_index_candidate_1, &check_result);						
							if (check_result == 0)
							{
								vertex_indices[counter] = vertex_index_candidate_1;
								latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
								longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
								++counter;
							}
						}
						if (counter != 5)
						{
							printf("Problem 15 in TRSK implementation detected.\n");
							exit(1);
						}
						for (int l = 0; l < 5; ++l)
						{
							if (vertex_indices[l] < 0 || vertex_indices[l] >= NUMBER_OF_DUAL_SCALARS_H)
							{
								printf("Problem 5 in TRSK implementation detected.\n");
								exit(1);
							}
							for (int m = l + 1; m < 5; ++m)
							{
								if (vertex_indices[l] == vertex_indices[m])
								{
									printf("Problem 23 in TRSK implementation detected.\n");
									exit(1);
								}
							}
						}
						retval = sort_edge_indices(latitude_vertices, longitude_vertices, 5, indices_resorted);
						for (int l = 0; l < 5; ++l)
						{
							for (int m = l + 1; m < 5; ++m)
							{
								if (indices_resorted[l] == indices_resorted[m] || indices_resorted[l] < 0 || indices_resorted[l] > 4)
								{
									printf("Problem 27 in TRSK implementation detected.\n");
									exit(1);
								}
							}
						}
						for (int l = 0; l < 5; ++l)
							vertex_indices_resorted[l] = vertex_indices[indices_resorted[l]];
						for (int l = 0; l < 5; ++l)
						{
							for (int m = 0; m < 5; ++m)
							{
								if ((from_index_dual[adjacent_vector_indices_h[6*to_index[i] + m]] == vertex_indices_resorted[l] && to_index_dual[adjacent_vector_indices_h[6*to_index[i] + m]] == vertex_indices_resorted[(l + 1)%5]) || (to_index_dual[adjacent_vector_indices_h[6*to_index[i] + m]] == vertex_indices_resorted[l] && from_index_dual[adjacent_vector_indices_h[6*to_index[i] + m]] == vertex_indices_resorted[(l + 1)%5]))
									edge_indices[l] = adjacent_vector_indices_h[6*to_index[i] + m];
							}
						}
						for (int l = 0; l < 5; ++l)
						{
							if (edge_indices[l] < 0 || edge_indices[l] >= NUMBER_OF_VECTORS_H)
							{
								printf("Problem 9 in TRSK implementation detected.\n");
								exit(1);
							}
							for (int m = l + 1; m < 5; ++m)
							{
								if (edge_indices[l] == edge_indices[m])
								{
									printf("Problem 20 in TRSK implementation detected.\n");
									exit(1);
								}
							}
						}
						for (int l = 0; l < 5; ++l)
						{
							latitude_edges[l] = latitude_vector[edge_indices[l]];
							longitude_edges[l] = longitude_vector[edge_indices[l]];
						}
						check_sum = 0;
						for (int l = 0; l < 5; ++l)	
						{
							if (l == 0)
								retval = calc_triangle_face(latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[4], longitude_edges[4], &triangle_0);
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
						for (int l = 0; l < 5; ++l)
						{
							if (edge_indices[l] == i)
								last_index = l;
							if (edge_indices[l] == trsk_modified_velocity_indices[10*i + k])
								first_index = (l + 1)%5;
								
						}
						if (k == 9)
							sum_of_weights = 0;
						else
							double_sum_gen(vector_of_areas, 5, first_index, last_index, &sum_of_weights);
						if (sum_of_weights < 0 || sum_of_weights/area[from_index[i]] > 1)
						{
							printf("Problem 36 in TRSK implementation detected\n");
							exit(1);
						}
					}
					else
					{
						int vertex_indices[6];
						int edge_indices[6];
						int indices_resorted[6];
						int vertex_indices_resorted[6];
						double latitude_vertices[6];
						double longitude_vertices[6];
						double latitude_edges[6];
						double longitude_edges[6];
						double vector_of_areas[6];
						for (int l = 0; l < 6; ++l)
							vertex_indices[l] = -1;
						counter = 0;
						for (int l = 0; l < 6; ++l)
						{
							vertex_index_candidate_0 = from_index_dual[adjacent_vector_indices_h[6*to_index[i] + l]];
							vertex_index_candidate_1 = to_index_dual[adjacent_vector_indices_h[6*to_index[i] + l]];
							retval = in_bool_calculator(vertex_indices, 6, vertex_index_candidate_0, &check_result);						
							if (check_result == 0)
							{
								vertex_indices[counter] = vertex_index_candidate_0;
								latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
								longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
								++counter;
							}
							retval = in_bool_calculator(vertex_indices, 6, vertex_index_candidate_1, &check_result);						
							if (check_result == 0)
							{
								vertex_indices[counter] = vertex_index_candidate_1;
								latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
								longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
								++counter;
							}
						}
						if (counter != 6)
						{
							printf("Problem 16 in TRSK implementation detected.\n");
							exit(1);
						}
						for (int l = 0; l < 6; ++l)
						{
							if (vertex_indices[l] < 0 || vertex_indices[l] >= NUMBER_OF_DUAL_SCALARS_H)
							{
								printf("Problem 6 in TRSK implementation detected.\n");
								exit(1);
							}
							for (int m = l + 1; m < 6; ++m)
							{
								if (vertex_indices[l] == vertex_indices[m])
								{
									printf("Problem 24 in TRSK implementation detected.\n");
									exit(1);
								}
							}
						}
						retval = sort_edge_indices(latitude_vertices, longitude_vertices, 6, indices_resorted);
						for (int l = 0; l < 6; ++l)
						{
							for (int m = l + 1; m < 6; ++m)
							{
								if (indices_resorted[l] == indices_resorted[m] || indices_resorted[l] < 0 || indices_resorted[l] > 5)
								{
									printf("Problem 28 in TRSK implementation detected.\n");
									exit(1);
								}
							}
						}
						for (int l = 0; l < 6; ++l)
							vertex_indices_resorted[l] = vertex_indices[indices_resorted[l]];
						for (int l = 0; l < 6; ++l)
						{
							for (int m = 0; m < 6; ++m)
							{
								if ((from_index_dual[adjacent_vector_indices_h[6*to_index[i] + m]] == vertex_indices_resorted[l] && to_index_dual[adjacent_vector_indices_h[6*to_index[i] + m]] == vertex_indices_resorted[(l + 1)%6]) || (to_index_dual[adjacent_vector_indices_h[6*to_index[i] + m]] == vertex_indices_resorted[l] && from_index_dual[adjacent_vector_indices_h[6*to_index[i] + m]] == vertex_indices_resorted[(l + 1)%6]))
									edge_indices[l] = adjacent_vector_indices_h[6*to_index[i] + m];
							}
						}
						for (int l = 0; l < 6; ++l)
						{
							if (edge_indices[l] < 0 || edge_indices[l] >= NUMBER_OF_VECTORS_H)
							{
								printf("Problem 10 in TRSK implementation detected.\n");
								exit(1);
							}
							for (int m = l + 1; m < 6; ++m)
							{
								if (edge_indices[l] == edge_indices[m])
								{
									printf("Problem 21 in TRSK implementation detected.\n");
									exit(1);
								}
							}
						}
						for (int l = 0; l < 6; ++l)
						{
							latitude_edges[l] = latitude_vector[edge_indices[l]];
							longitude_edges[l] = longitude_vector[edge_indices[l]];
						}
						check_sum = 0;
						for (int l = 0; l < 6; ++l)	
						{
							if (l == 0)
								retval = calc_triangle_face(latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[5], longitude_edges[5], &triangle_0);
							else
								retval = calc_triangle_face(latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[l - 1], longitude_edges[l - 1], &triangle_0);
							retval = calc_triangle_face(latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[l], longitude_edges[l], &triangle_1);
							vector_of_areas[l] = pow(RADIUS + z_scalar[to_index[i]], 2)*(triangle_0 + triangle_1);
							check_sum += vector_of_areas[l];
						}
						if (fabs(check_sum/area[to_index[i]] - 1) > 0.001)
						{
							printf("Problem 33 in TRSK implementation detected.\n");
							exit(1);
						}
						for (int l = 0; l < 6; ++l)
						{
							if (edge_indices[l] == i)
								last_index = l;
							if (edge_indices[l] == trsk_modified_velocity_indices[10*i + k])
								first_index = (l + 1)%6;
						}
						double_sum_gen(vector_of_areas, 6, first_index, last_index, &sum_of_weights);
						if (sum_of_weights < 0 || sum_of_weights/area[from_index[i]] > 1)
						{
							printf("Problem 37 in TRSK implementation detected\n");
							exit(1);
						}
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
	// free(adjacent_scalar_indices_for_cross); uncomment for old Coriolis reconstruction
    double weight_prefactor;
    for (int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
    {
        weight_prefactor = 2.0/6.0;
        if (i < NUMBER_OF_PENTAGONS)
            weight_prefactor = 2.0/5.0;
        for (int j = 0; j < 6; ++j)
        {
            recov_ver_index[6*i + j] = adjacent_vector_indices_h[6*i + j];
            recov_ver_0_pri_weight[6*i + j] = weight_prefactor*cos(direction[recov_ver_index[6*i + j]]);
            recov_ver_0_curl_weight[6*i + j] = weight_prefactor*cos(direction_dual[recov_ver_index[6*i + j]]);
            recov_ver_1_pri_weight[6*i + j] = weight_prefactor*sin(direction[recov_ver_index[6*i + j]]);
            recov_ver_1_curl_weight[6*i + j] = weight_prefactor*sin(direction_dual[recov_ver_index[6*i + j]]);
        }
        if (i < NUMBER_OF_PENTAGONS)
        {
            recov_ver_index[6*i + 5] = 0;
            recov_ver_0_pri_weight[6*i + 5] = 0;
            recov_ver_1_pri_weight[6*i + 5] = 0;
            recov_ver_0_curl_weight[6*i + 5] = 0;
            recov_ver_1_curl_weight[6*i + 5] = 0;
        }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_H; ++i)
    {
        sign = 1;
        find_angle_change(direction_dual[i], direction[i], &direction_change);
        if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
            sign = -1;
        h_curl_indices[4*i + 0] = i + NUMBER_OF_VECTORS_PER_LAYER;
        h_curl_signs[4*i + 0] = sign;
        if (sign == 1)
            h_curl_indices[4*i + 1] = to_index[i];
        else
            h_curl_indices[4*i + 1] = from_index[i];
        h_curl_signs[4*i + 1] = 1;
        h_curl_indices[4*i + 2] = i;
        h_curl_signs[4*i + 2] = -sign;
        if (sign == 1)
            h_curl_indices[4*i + 3] = from_index[i];
        else
            h_curl_indices[4*i + 3] = to_index[i];
        h_curl_signs[4*i + 3] = -1;
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
    {
    	counter = 0;
    	for (int j = 0; j < NUMBER_OF_DUAL_VECTORS_H; ++j)
    	{
    		if (h_curl_indices[4*j + 1] == i || h_curl_indices[4*j + 3] == i)
    			++counter;
    	}
    	if (i < NUMBER_OF_PENTAGONS && counter != 5)
    	{
    		printf("Error in h_curl_indices, position 0.\n");
    		exit(1);
		}
    	if (i >= NUMBER_OF_PENTAGONS && counter != 6)
    	{
    		printf("Error in h_curl_indices, position 1.\n");
    		exit(1);
		}
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
    	counter = 0;
    	for (int j = 0; j < NUMBER_OF_DUAL_VECTORS_H; ++j)
    	{
    		if (h_curl_indices[4*j + 0] == i + NUMBER_OF_VECTORS_PER_LAYER || h_curl_indices[4*j + 2] == i)
    		{
    			++counter;
    			if (h_curl_indices[4*j + 0] == i + NUMBER_OF_VECTORS_PER_LAYER && h_curl_indices[4*j + 2] == i)
    				++counter;
			}
    	}
    	if (counter != 2)
    	{
    		printf("Error in h_curl_indices, position 2.\n");
    		exit(1);
		}
    }
    free(direction_dual);
    double area_0, area_1, area_ratio;
    for (int i = 0; i < NUMBER_OF_DUAL_H_VECTORS + NUMBER_OF_H_VECTORS; ++i)
    {
    	layer_index = i/(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H);
    	h_index = i - layer_index*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H);
    	if (h_index < NUMBER_OF_DUAL_VECTORS_H)
    		area_dual[i] = area_dual_pre[layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + h_index];
    	else
    	{
    		primal_vector_index = NUMBER_OF_VECTORS_V + h_index - NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_VECTORS_PER_LAYER;
    		dual_vector_index = NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + to_index_dual[h_index - NUMBER_OF_DUAL_VECTORS_H];
    		area_0 = area_dual_pre[dual_vector_index];
    		area_rescale_factor = pow((RADIUS + z_vector[primal_vector_index])/(RADIUS + z_vector_dual[dual_vector_index]), 2);
    		area_0 = area_rescale_factor*area_0;
    		dual_vector_index = NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER + from_index_dual[h_index - NUMBER_OF_DUAL_VECTORS_H];
    		area_1 = area_dual_pre[dual_vector_index];
    		area_rescale_factor = pow((RADIUS + z_vector[primal_vector_index])/(RADIUS + z_vector_dual[dual_vector_index]), 2);
    		area_1 = area_rescale_factor*area_1;
    		area_ratio = area_0/area_1;
    		if (fabs(area_ratio - 1) > 0.3)
    		{
    			printf("Unrealistic area_ratio in rhombus area calculation, position 0.\n");
    			exit(1);
			}
    		area_dual[i] = area_0 + area_1;
    	}
		if (area_dual[i] <= 0)
		{
			printf("It is area_dual <= 0 at some point.\n");
			exit(1);
		}
    }
    double check_area, area_sum, mean_z;
    for (int i = 0; i < NUMBER_OF_LAYERS; ++i)
    {
    	mean_z = 0;
    	for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
    		mean_z += z_vector[NUMBER_OF_VECTORS_V + i*NUMBER_OF_VECTORS_PER_LAYER + j];
    	mean_z = mean_z/NUMBER_OF_VECTORS_H;
    	area_sum = 0;
    	for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
    		area_sum += 1.0/3*area_dual[NUMBER_OF_DUAL_VECTORS_H + i*(NUMBER_OF_VECTORS_H + NUMBER_OF_DUAL_VECTORS_H) + j];
    	check_area = 4*M_PI*pow(RADIUS + mean_z, 2);
    	area_ratio = check_area/area_sum;
    	if (fabs(area_ratio - 1) > 0.00001)
    	{
    		printf("Unrealistic area_ratio in rhombus area calculation, position 1.\n");
    		exit(1);
    	}
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_H; ++i)
    {
    	area_ratio = area_dual[i]/area_dual[i + NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H];
    	if (fabs(area_ratio - 0.5) > 0.0001)
    	{
    		printf("Unrealistic value in area_ratio of area_dual, position 0.\n");
    		exit(1);
    	}
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_H; ++i)
    {
    	area_ratio = area_dual[i + (NUMBER_OF_LAYERS - 1)*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H)]/area_dual[i + (NUMBER_OF_LAYERS - 1)*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H) + NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H];
    	if (fabs(area_ratio - 2) > 0.001)
    	{
    		printf("Unrealistic value in area_ratio of area_dual, position 1.\n");
    		exit(1);
    	}
    }
    int indices_list_pre[6];
    int signs_list_pre[6];
    int indices_list[4];
    int signs_list[4];
    int double_indices[2];
    double_indices[0] = -1;
    double_indices[1] = -1;
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
    	for (int j = 0; j < 3; ++j)
    	{
			indices_list_pre[j] = vorticity_indices_pre[3*to_index_dual[i] + j];
			signs_list_pre[j] = vorticity_signs_pre[3*to_index_dual[i] + j];
    	}
    	for (int j = 0; j < 3; ++j)
    	{
			indices_list_pre[3 + j] = vorticity_indices_pre[3*from_index_dual[i] + j];
			signs_list_pre[3 + j] = vorticity_signs_pre[3*from_index_dual[i] + j];
    	}
		for (int j = 0; j < 6; ++j)
		{
			for (int k = j + 1; k < 6; ++k)
			{
				if (indices_list_pre[j] == indices_list_pre[k])
				{
					double_indices[0] = j;
					double_indices[1] = k;
				}
			}
		}
		counter = 0;
    	for (int j = 0; j < 6; ++j)
    	{
    		if (j != double_indices[0] && j != double_indices[1])
    		{
				indices_list[counter] = indices_list_pre[j];
				signs_list[counter] = signs_list_pre[j];
				counter++;
			}
		}
		if (counter != 4)
		{
			printf("Error in vorticity_indices and vortictiy_signs creation from vorticity_indices_pre and vortictiy_signs_pre, position 1.\n");
			exit(1);
		}
    	for (int j = 0; j < 4; ++j)
    	{
			vorticity_indices[4*i + j] = indices_list[j];
			vorticity_signs[4*i + j] = signs_list[j];
			if (vorticity_signs[4*i + j] != 1 && vorticity_signs[4*i + j] != -1)
			{
				printf("Error in vorticity_indices and vortictiy_signs creation from vorticity_indices_pre and vortictiy_signs_pre, position 2.");
				exit(1);
			}
			if (vorticity_indices[4*i + j] >= NUMBER_OF_VECTORS_H || vorticity_indices[4*i + j] < 0)
			{
				printf("Error in vorticity_indices and vortictiy_signs creation from vorticity_indices_pre and vortictiy_signs_pre, position 3.");
				exit(1);
			}
		}
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
    	counter = 0;
    	for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
    	{
    		for (int k = 0; k < 4; ++k)
    		{
    			if (vorticity_indices[4*j + k] == i)
    				++counter;
    		}
    	}
    	if (counter != 4)
    	{
    		printf("Error in vorticity_indices, position 0.\n");
    		exit(1);
    	}
    }
    int ncid_g_prop;
    if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid_g_prop)))
        ERR(retval);
    free(OUTPUT_FILE);
    int latitude_scalar_id, longitude_scalar_id, direction_id, latitude_vector_id, longitude_vector_id, latitude_scalar_dual_id, longitude_scalar_dual_id, z_scalar_id, z_vector_id, normal_distance_id, gravity_id, volume_id, area_id, recov_hor_par_curl_weight_id, recov_hor_ver_curl_weight_id, trsk_modified_weights_id, recov_ver_0_pri_weight_id, recov_ver_0_curl_weight_id, recov_ver_1_pri_weight_id, recov_ver_1_curl_weight_id, z_vector_dual_id, normal_distance_dual_id, area_dual_id, f_vec_id, to_index_id, from_index_id, adjacent_vector_indices_h_id, vorticity_indices_id, h_curl_indices_id, recov_hor_par_curl_index_id, recov_hor_ver_curl_index_id, trsk_modified_velocity_indices_id, trsk_modified_curl_indices_id, recov_hor_ver_pri_index_id, recov_ver_index_id, adjacent_signs_h_id, vorticity_signs_id, h_curl_signs_id, vector_curl_one_layer_dimid, scalar_dimid, scalar_h_dimid, scalar_dual_h_dimid, vector_dimid, scalar_h_dimid_6, vector_h_dimid, vector_h_dimid_11, vector_h_dimid_10, vector_h_dimid_2, vector_h_dimid_4, vector_v_dimid_6, vector_dual_dimid, vector_dual_h_dimid, vector_dual_v_dimid_3, vector_dual_h_dimid_4, adjacent_scalar_indices_dual_h_id, gravity_potential_id, vertical_contravar_unit_id, vertical_contravar_unit_dimid, adjacent_vector_indices_dual_h_id, scalar_dual_h_dimid_3, vector_dual_area_dimid, e_kin_weights_id, scalar_14_dimid, e_kin_indices_id;
    printf("Starting to write to output file ... ");
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_index", NUMBER_OF_SCALARS, &scalar_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_14_index", 14*NUMBER_OF_SCALARS, &scalar_14_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_h_index", NUMBER_OF_SCALARS_H, &scalar_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_dual_h_index", NUMBER_OF_DUAL_SCALARS_H, &scalar_dual_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_dual_h_3_index", 3*NUMBER_OF_DUAL_SCALARS_H, &scalar_dual_h_dimid_3)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index", NUMBER_OF_VECTORS, &vector_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_index", NUMBER_OF_VECTORS_H, &vector_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_h_6_index", 6*NUMBER_OF_SCALARS_H, &scalar_h_dimid_6)))
        ERR(retval);
	if ((retval = nc_def_dim(ncid_g_prop, "vector_h_10_index", 10*NUMBER_OF_VECTORS_H, &vector_h_dimid_10)))
	    ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_11_index", 11*NUMBER_OF_VECTORS_H, &vector_h_dimid_11)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_2_index", 2*NUMBER_OF_VECTORS_H, &vector_h_dimid_2)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_4_index", 4*NUMBER_OF_VECTORS_H, &vector_h_dimid_4)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_v_6_index", 6*NUMBER_OF_VECTORS_V, &vector_v_dimid_6)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_contravar_unit_index", 3*NUMBER_OF_VECTORS_V*(NUMBER_OF_ORO_LAYERS + 1), &vertical_contravar_unit_dimid)))
        ERR(retval);        
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_h_dual", NUMBER_OF_DUAL_VECTORS_H, &vector_dual_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_curl_one_layer_dimid", NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_VECTORS_H, &vector_curl_one_layer_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_dual", NUMBER_OF_DUAL_VECTORS, &vector_dual_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_dual_area", NUMBER_OF_DUAL_H_VECTORS + NUMBER_OF_H_VECTORS, &vector_dual_area_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_dual_v_3_index", 3*NUMBER_OF_DUAL_VECTORS_V, &vector_dual_v_dimid_3)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_dual_h_4_index", 4*NUMBER_OF_DUAL_VECTORS_H, &vector_dual_h_dimid_4)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "latitude_scalar", NC_DOUBLE, 1, &scalar_h_dimid, &latitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "longitude_scalar", NC_DOUBLE, 1, &scalar_h_dimid, &longitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "latitude_scalar_dual", NC_DOUBLE, 1, &scalar_dual_h_dimid, &latitude_scalar_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "longitude_scalar_dual", NC_DOUBLE, 1, &scalar_dual_h_dimid, &longitude_scalar_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_scalar", NC_DOUBLE, 1, &scalar_dimid, &z_scalar_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_scalar_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "gravity_potential", NC_DOUBLE, 1, &scalar_dimid, &gravity_potential_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, gravity_potential_id, "units", strlen("m^2/s^2"), "m^2/s^2")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vertical_contravar_unit", NC_DOUBLE, 1, &vertical_contravar_unit_dimid, &vertical_contravar_unit_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_surface", NC_DOUBLE, 1, &scalar_h_dimid, &z_surface_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_surface_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_vector", NC_DOUBLE, 1, &vector_dimid, &z_vector_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_vector_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "normal_distance", NC_DOUBLE, 1, &vector_dimid, &normal_distance_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, normal_distance_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "gravity", NC_DOUBLE, 1, &vector_dimid, &gravity_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, gravity_id, "units", strlen("m/s^2"), "m/s^2")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "volume", NC_DOUBLE, 1, &scalar_dimid, &volume_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, volume_id, "units", strlen("m^3"), "m^3")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "area", NC_DOUBLE, 1, &vector_dimid, &area_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, area_id, "units", strlen("m^2"), "m^2")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_par_curl_weight", NC_DOUBLE, 1, &vector_h_dimid_2, &recov_hor_par_curl_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_ver_curl_weight", NC_DOUBLE, 1, &vector_h_dimid_2, &recov_hor_ver_curl_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "trsk_modified_weights", NC_DOUBLE, 1, &vector_h_dimid_10, &trsk_modified_weights_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_0_curl_weight", NC_DOUBLE, 1, &vector_v_dimid_6, &recov_ver_0_curl_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_0_pri_weight", NC_DOUBLE, 1, &vector_v_dimid_6, &recov_ver_0_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_1_pri_weight", NC_DOUBLE, 1, &vector_v_dimid_6, &recov_ver_1_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_1_curl_weight", NC_DOUBLE, 1, &vector_v_dimid_6, &recov_ver_1_curl_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_vector_dual", NC_DOUBLE, 1, &vector_dual_dimid, &z_vector_dual_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_vector_dual_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "normal_distance_dual", NC_DOUBLE, 1, &vector_dual_dimid, &normal_distance_dual_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, normal_distance_dual_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "area_dual", NC_DOUBLE, 1, &vector_dual_area_dimid, &area_dual_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, area_dual_id, "units", strlen("m^2"), "m^2")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "f_vec", NC_DOUBLE, 1, &vector_curl_one_layer_dimid, &f_vec_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, f_vec_id, "units", strlen("1/s"), "1/s")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "direction", NC_DOUBLE, 1, &vector_h_dimid, &direction_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "latitude_vector", NC_DOUBLE, 1, &vector_h_dimid, &latitude_vector_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "longitude_vector", NC_DOUBLE, 1, &vector_h_dimid, &longitude_vector_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "e_kin_weights", NC_DOUBLE, 1, &scalar_14_dimid, &e_kin_weights_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "e_kin_indices", NC_INT, 1, &scalar_14_dimid, &e_kin_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "to_index", NC_INT, 1, &vector_h_dimid, &to_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "from_index", NC_INT, 1, &vector_h_dimid, &from_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_vector_indices_h", NC_INT, 1, &scalar_h_dimid_6, &adjacent_vector_indices_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_vector_indices_dual_h", NC_INT, 1, &scalar_dual_h_dimid_3, &adjacent_vector_indices_dual_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_indices", NC_INT, 1, &vector_h_dimid_4, &vorticity_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "h_curl_indices", NC_INT, 1, &vector_dual_h_dimid_4, &h_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_par_curl_index", NC_INT, 1, &vector_h_dimid_2, &recov_hor_par_curl_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_ver_curl_index", NC_INT, 1, &vector_h_dimid_2, &recov_hor_ver_curl_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "trsk_modified_velocity_indices", NC_INT, 1, &vector_h_dimid_10, &trsk_modified_velocity_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "trsk_modified_curl_indices", NC_INT, 1, &vector_h_dimid_10, &trsk_modified_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_ver_pri_index", NC_INT, 1, &vector_h_dimid_4, &recov_hor_ver_pri_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_index", NC_INT, 1, &vector_v_dimid_6, &recov_ver_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_scalar_indices_dual_h", NC_INT, 1, &vector_dual_v_dimid_3, &adjacent_scalar_indices_dual_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_signs_h", NC_INT, 1, &scalar_h_dimid_6, &adjacent_signs_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_signs", NC_INT, 1, &vector_h_dimid_4, &vorticity_signs_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "h_curl_signs", NC_INT, 1, &vector_dual_h_dimid_4, &h_curl_signs_id)))
        ERR(retval);
    if ((retval = nc_enddef(ncid_g_prop)))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, latitude_scalar_id, &latitude_scalar[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, longitude_scalar_id, &longitude_scalar[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, latitude_scalar_dual_id, &latitude_scalar_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, longitude_scalar_dual_id, &longitude_scalar_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, z_scalar_id, &z_scalar[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, gravity_potential_id, &gravity_potential[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, vertical_contravar_unit_id, &vertical_contravar_unit[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, z_surface_id, &z_surface[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, z_vector_id, &z_vector[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, normal_distance_id, &normal_distance[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, gravity_id, &gravity[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, volume_id, &volume[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, area_id, &area[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, e_kin_weights_id, &e_kin_weights[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_hor_par_curl_weight_id, &recov_hor_par_curl_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_hor_ver_curl_weight_id, &recov_hor_ver_curl_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, trsk_modified_weights_id, &trsk_modified_weights[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_0_pri_weight_id, &recov_ver_0_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_0_curl_weight_id, &recov_ver_0_curl_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_1_pri_weight_id, &recov_ver_1_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_1_curl_weight_id, &recov_ver_1_curl_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, z_vector_dual_id, &z_vector_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, normal_distance_dual_id, &normal_distance_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, area_dual_id, &area_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, f_vec_id, &f_vec[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, direction_id, &direction[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, latitude_vector_id, &latitude_vector[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, longitude_vector_id, &longitude_vector[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, to_index_id, &to_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, from_index_id, &from_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, adjacent_vector_indices_h_id, &adjacent_vector_indices_h[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, adjacent_vector_indices_dual_h_id, &adjacent_vector_indices_dual_h[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, e_kin_indices_id, &e_kin_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, vorticity_indices_id, &vorticity_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, h_curl_indices_id, &h_curl_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, recov_hor_par_curl_index_id, &recov_hor_par_curl_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, recov_hor_ver_curl_index_id, &recov_hor_ver_curl_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, trsk_modified_velocity_indices_id, &trsk_modified_velocity_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, trsk_modified_curl_indices_id, &trsk_modified_curl_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, recov_hor_ver_pri_index_id, &recov_hor_ver_pri_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, recov_ver_index_id, &recov_ver_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, adjacent_scalar_indices_dual_h_id, &adjacent_scalar_indices_dual_h[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, adjacent_signs_h_id, &adjacent_signs_h[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, vorticity_signs_id, &vorticity_signs[0])))
        ERR(retval);
    if ((retval = nc_put_var_int(ncid_g_prop, h_curl_signs_id, &h_curl_signs[0])))
        ERR(retval);
    if ((retval = nc_close(ncid_g_prop)))
        ERR(retval);
    printf("finished.\n");
    free(e_kin_weights);
    free(e_kin_indices);
    free(adjacent_vector_indices_dual_h);
    free(vertical_contravar_unit);
    free(gravity_potential);
    free(adjacent_scalar_indices_dual_h);
    free(latitude_vector);
    free(longitude_vector);
    free(direction);
    free(latitude_scalar);
    free(longitude_scalar);
    free(z_scalar);
    free(z_vector);
    free(normal_distance);
    free(gravity);
    free(volume);
    free(area);
    free(recov_hor_par_curl_weight);
    free(recov_hor_ver_curl_weight);
    free(trsk_modified_weights);
    free(recov_ver_0_pri_weight);
    free(recov_ver_0_curl_weight);
    free(recov_ver_1_pri_weight);
    free(recov_ver_1_curl_weight);
    free(latitude_scalar_dual);
    free(longitude_scalar_dual);
    free(z_scalar_dual);
    free(latitude_vector_dual);
    free(z_vector_dual);
    free(normal_distance_dual);
    free(f_vec);
    free(to_index);
    free(from_index);
    free(to_index_dual);
    free(from_index_dual);
    free(adjacent_vector_indices_h);
    free(vorticity_indices_pre);
    free(vorticity_indices);
    free(h_curl_indices);
    free(recov_hor_par_curl_index);
    free(recov_hor_ver_curl_index);
    free(trsk_modified_velocity_indices);
    free(trsk_modified_curl_indices);
    free(recov_hor_ver_pri_index);
    free(recov_ver_index);
    free(adjacent_signs_h);
    free(vorticity_signs_pre);
    free(vorticity_signs);
    free(h_curl_signs);
    free(area_dual_pre);
    free(area_dual);
	free(z_surface);
    return 0;
}






