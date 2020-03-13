#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include "/home/max/my_code/game/core/src/enum_and_typedefs.h"
#include "/home/max/my_code/c_source_libs/geos/geos.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define FILE_NAME "nc_files/res_3_oro_0_geo_prop.nc"

int main(int argc, char *argv[])
{
    const double ATMOS_HEIGHT = SCALE_HEIGHT*log(NUMBER_OF_LEVELS);
    double *latitude_ico = malloc(12*sizeof(double));
    latitude_ico[0] = M_PI/2;
    latitude_ico[1] = M_PI/6;
    latitude_ico[2] = M_PI/6;
    latitude_ico[3] = M_PI/6;
    latitude_ico[4] = M_PI/6;
    latitude_ico[5] = M_PI/6;
    latitude_ico[6] = -M_PI/6;
    latitude_ico[7] = -M_PI/6;
    latitude_ico[8] = -M_PI/6;
    latitude_ico[9] = -M_PI/6;
    latitude_ico[10] = -M_PI/6;
    latitude_ico[11] = -M_PI/2;
    double *longitude_ico = malloc(12*sizeof(double));
    longitude_ico[0] = 0;
    longitude_ico[1] = 0;
    longitude_ico[2] = 1*2*M_PI/5;
    longitude_ico[3] = 2*2*M_PI/5;
    longitude_ico[4] = 3*2*M_PI/5;
    longitude_ico[5] = 4*2*M_PI/5;
    longitude_ico[6] = 2*M_PI/10;
    longitude_ico[7] = 2*M_PI/10 + 1*2*M_PI/5;
    longitude_ico[8] = 2*M_PI/10 + 2*2*M_PI/5;
    longitude_ico[9] = 2*M_PI/10 + 3*2*M_PI/5;
    longitude_ico[10] = 2*M_PI/10 + 4*2*M_PI/5;    
    longitude_ico[11] = 0;
    int edge_vertices[NUMBER_OF_EDGES][2];
    edge_vertices[0][0] = 0;
    edge_vertices[0][1] = 1;
    edge_vertices[1][0] = 0;
    edge_vertices[1][1] = 2;
    edge_vertices[2][0] = 0;
    edge_vertices[2][1] = 3;
    edge_vertices[3][0] = 0;
    edge_vertices[3][1] = 4;
    edge_vertices[4][0] = 0;
    edge_vertices[4][1] = 5;
    edge_vertices[5][0] = 1;
    edge_vertices[5][1] = 2;
    edge_vertices[6][0] = 2;
    edge_vertices[6][1] = 3;
    edge_vertices[7][0] = 3;
    edge_vertices[7][1] = 4;
    edge_vertices[8][0] = 4;
    edge_vertices[8][1] = 5;
    edge_vertices[9][0] = 5;
    edge_vertices[9][1] = 1;
    edge_vertices[10][0] = 1;
    edge_vertices[10][1] = 6;
    edge_vertices[11][0] = 2;
    edge_vertices[11][1] = 6;
    edge_vertices[12][0] = 2;
    edge_vertices[12][1] = 7;
    edge_vertices[13][0] = 3;
    edge_vertices[13][1] = 7;
    edge_vertices[14][0] = 3;
    edge_vertices[14][1] = 8;
    edge_vertices[15][0] = 4;
    edge_vertices[15][1] = 8;
    edge_vertices[16][0] = 4;
    edge_vertices[16][1] = 9;
    edge_vertices[17][0] = 5;
    edge_vertices[17][1] = 9;
    edge_vertices[18][0] = 5;
    edge_vertices[18][1] = 10;
    edge_vertices[19][0] = 1;
    edge_vertices[19][1] = 10;
    edge_vertices[20][0] = 10;
    edge_vertices[20][1] = 6;
    edge_vertices[21][0] = 6;
    edge_vertices[21][1] = 7;
    edge_vertices[22][0] = 7;
    edge_vertices[22][1] = 8;
    edge_vertices[23][0] = 8;
    edge_vertices[23][1] = 9;
    edge_vertices[24][0] = 9;
    edge_vertices[24][1] = 10;
    edge_vertices[25][0] = 6;
    edge_vertices[25][1] = 11;
    edge_vertices[26][0] = 7;
    edge_vertices[26][1] = 11;
    edge_vertices[27][0] = 8;
    edge_vertices[27][1] = 11;
    edge_vertices[28][0] = 9;
    edge_vertices[28][1] = 11;
    edge_vertices[29][0] = 10;
    edge_vertices[29][1] = 11;
    int face_vertices[20][3];
    face_vertices[0][0] = 0;
    face_vertices[0][1] = 1;
    face_vertices[0][2] = 2;
    face_vertices[1][0] = 0;
    face_vertices[1][1] = 2;
    face_vertices[1][2] = 3;
    face_vertices[2][0] = 0;
    face_vertices[2][1] = 3;
    face_vertices[2][2] = 4;
    face_vertices[3][0] = 0;
    face_vertices[3][1] = 4;
    face_vertices[3][2] = 5;
    face_vertices[4][0] = 0;
    face_vertices[4][1] = 5;
    face_vertices[4][2] = 1;
    face_vertices[5][0] = 1;
    face_vertices[5][1] = 10;
    face_vertices[5][2] = 6;
    face_vertices[6][0] = 6;
    face_vertices[6][1] = 2;
    face_vertices[6][2] = 1;
    face_vertices[7][0] = 2;
    face_vertices[7][1] = 6;
    face_vertices[7][2] = 7;
    face_vertices[8][0] = 7;
    face_vertices[8][1] = 3;
    face_vertices[8][2] = 2;
    face_vertices[9][0] = 3;
    face_vertices[9][1] = 7;
    face_vertices[9][2] = 8;
    face_vertices[10][0] = 8;
    face_vertices[10][1] = 4;
    face_vertices[10][2] = 3;
    face_vertices[11][0] = 4;
    face_vertices[11][1] = 8;
    face_vertices[11][2] = 9;
    face_vertices[12][0] = 9;
    face_vertices[12][1] = 5;
    face_vertices[12][2] = 4;
    face_vertices[13][0] = 5;
    face_vertices[13][1] = 9;
    face_vertices[13][2] = 10;
    face_vertices[14][0] = 10;
    face_vertices[14][1] = 1;
    face_vertices[14][2] = 5;
    face_vertices[15][0] = 11;
    face_vertices[15][1] = 6;
    face_vertices[15][2] = 10;
    face_vertices[16][0] = 11;
    face_vertices[16][1] = 7;
    face_vertices[16][2] = 6;
    face_vertices[17][0] = 11;
    face_vertices[17][1] = 8;
    face_vertices[17][2] = 7;
    face_vertices[18][0] = 11;
    face_vertices[18][1] = 9;
    face_vertices[18][2] = 8;
    face_vertices[19][0] = 11;
    face_vertices[19][1] = 10;
    face_vertices[19][2] = 9;
    int face_edges[20][3];
    face_edges[0][0] = 0;
    face_edges[0][1] = 5;
    face_edges[0][2] = 1;
    face_edges[1][0] = 1;
    face_edges[1][1] = 6;
    face_edges[1][2] = 2;
    face_edges[2][0] = 2;
    face_edges[2][1] = 7;
    face_edges[2][2] = 3;
    face_edges[3][0] = 3;
    face_edges[3][1] = 8;
    face_edges[3][2] = 4;
    face_edges[4][0] = 4;
    face_edges[4][1] = 9;
    face_edges[4][2] = 0;
    face_edges[5][0] = 19;
    face_edges[5][1] = 20;
    face_edges[5][2] = 10;
    face_edges[6][0] = 11;
    face_edges[6][1] = 5;
    face_edges[6][2] = 10;
    face_edges[7][0] = 11;
    face_edges[7][1] = 21;
    face_edges[7][2] = 12;
    face_edges[8][0] = 13;
    face_edges[8][1] = 6;
    face_edges[8][2] = 12;
    face_edges[9][0] = 13;
    face_edges[9][1] = 22;
    face_edges[9][2] = 14;
    face_edges[10][0] = 15;
    face_edges[10][1] = 7;
    face_edges[10][2] = 14;
    face_edges[11][0] = 15;
    face_edges[11][1] = 23;
    face_edges[11][2] = 16;
    face_edges[12][0] = 17;
    face_edges[12][1] = 8;
    face_edges[12][2] = 16;
    face_edges[13][0] = 17;
    face_edges[13][1] = 24;
    face_edges[13][2] = 18;
    face_edges[14][0] = 19;
    face_edges[14][1] = 9;
    face_edges[14][2] = 18;
    face_edges[15][0] = 25;
    face_edges[15][1] = 20;
    face_edges[15][2] = 29;
    face_edges[16][0] = 26;
    face_edges[16][1] = 21;
    face_edges[16][2] = 25;
    face_edges[17][0] = 27;
    face_edges[17][1] = 22;
    face_edges[17][2] = 26;
    face_edges[18][0] = 28;
    face_edges[18][1] = 23;
    face_edges[18][2] = 27;
    face_edges[19][0] = 29;
    face_edges[19][1] = 24;
    face_edges[19][2] = 28;
    double *latitude_scalar = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *z_scalar = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *z_vector = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *normal_distance = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *gravity = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *volume = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *area = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *recov_hor_par_dual_weight = malloc(11*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_hor_ver_dual_weight = malloc(2*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_hor_par_pri_weight = malloc(2*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_hor_ver_pri_weight = malloc(4*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_ver_1_pri_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *recov_ver_1_dual_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *recov_ver_2_pri_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *recov_ver_2_dual_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *latitude_scalar_dual = malloc(NUMBER_OF_DUAL_SCALARS_H*sizeof(double));
    double *longitude_scalar_dual = malloc(NUMBER_OF_DUAL_SCALARS_H*sizeof(double));
    double *z_scalar_dual = malloc(NUMBER_OF_DUAL_SCALARS*sizeof(double));
    double *latitude_vector_dual = malloc(NUMBER_OF_DUAL_VECTORS_PER_LAYER*sizeof(double));
    double *z_vector_dual = malloc(NUMBER_OF_DUAL_VECTORS*sizeof(double));
    double *normal_distance_dual = malloc(NUMBER_OF_DUAL_VECTORS*sizeof(double));
    double *area_dual = malloc(NUMBER_OF_DUAL_VECTORS*sizeof(double));
    double *f_vec = malloc(NUMBER_OF_DUAL_VECTORS_PER_LAYER*sizeof(double));
    long *to_index = malloc(NUMBER_OF_VECTORS_H*sizeof(long));
    long *from_index = malloc(NUMBER_OF_VECTORS_H*sizeof(long));
    long *recov_ver_2_dual_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(long));
    long *recov_ver_1_dual_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(long));
    long *recov_ver_2_pri_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(long));
    long *recov_hor_ver_pri_index = malloc(4*NUMBER_OF_VECTORS_H*sizeof(long));
    long *recov_ver_1_pri_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(long));
    long *recov_hor_ver_dual_index = malloc(2*NUMBER_OF_VECTORS_H*sizeof(long));
    long *recov_hor_par_pri_index = malloc(2*NUMBER_OF_VECTORS_H*sizeof(long));
    long *recov_hor_par_dual_index = malloc(11*NUMBER_OF_VECTORS_H*sizeof(long));
    long *adjacent_vector_indices_h = malloc(6*NUMBER_OF_SCALARS_H*sizeof(long));
    long *vorticity_indices = malloc(6*NUMBER_OF_SCALARS_H*sizeof(long));
    long *h_curl_indices = malloc(4*NUMBER_OF_VECTORS_H*sizeof(long));
    long *to_index_dual = malloc(NUMBER_OF_DUAL_VECTORS_H*sizeof(long));
    long *from_index_dual = malloc(NUMBER_OF_DUAL_VECTORS_H*sizeof(long));
    long *vorticity_indices_dual = malloc(3*NUMBER_OF_DUAL_VECTORS_V*sizeof(long));
    long *h_curl_indices_dual = malloc(4*NUMBER_OF_DUAL_VECTORS_H*sizeof(long));
    short *adjacent_signs_h = malloc(6*NUMBER_OF_SCALARS_H*sizeof(short));
    short *vorticity_signs = malloc(6*NUMBER_OF_SCALARS_H*sizeof(short));
    short *h_curl_signs = malloc(4*NUMBER_OF_VECTORS_H*sizeof(short));
    short *vector_product_sign = malloc(NUMBER_OF_VECTORS_H*sizeof(short));
    short *vorticity_signs_dual = malloc(3*NUMBER_OF_DUAL_VECTORS_V*sizeof(short));
    short *h_curl_signs_dual = malloc(4*NUMBER_OF_DUAL_VECTORS_H*sizeof(short));
    double lat_edge_1, lon_edge_1, lat_edge_2, lon_edge_2, point_frac, lat_1, lon_1, lat_2, lon_2, rel_on_line, lat_res, lon_res, base_area, base_distance, radius_1, radius_2, z_1, z_2, r_1, r_2, parallel_distance;
    int small_triangle_edge_index, reverse_0, reverse_1, reverse_2, face_index;
    long h_index, edge_index, on_face_index, layer_index, level_index, inner_index, upper_index, lower_index, check, edge_0, edge_1, edge_2, edge_index_1, edge_index_2, coord_1_points_amount, index_check, j, coord_1, coord_2, vertex_index_0, vertex_index_1, vertex_index_2, index_1, index_2, index_3, points_right, triangle_on_face_index, floor_index, on_edge_index, point_0, point_1, point_2, point_3, dual_scalar_index, counter, primal_vector_index, dual_vector_index;
    for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
    {
        if (i < NUMBER_OF_PENTAGONS)
        {
            latitude_scalar[i] = latitude_ico[i];
            longitude_scalar[i] = longitude_ico[i];
        }
        else if (i < NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES)
        {
            edge_index = (i - NUMBER_OF_PENTAGONS)/POINTS_PER_EDGE;
            point_frac = (1.0 + i - (NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE))/(1.0 + POINTS_PER_EDGE);
            vertex_index_0 = edge_vertices[edge_index][0];
            vertex_index_1 = edge_vertices[edge_index][1];
            find_geodetic(latitude_ico[vertex_index_0], longitude_ico[vertex_index_0], latitude_ico[vertex_index_1], longitude_ico[vertex_index_1], point_frac, &lat_res, &lon_res);
            latitude_scalar[i] = lat_res;
            longitude_scalar[i] = lon_res;
        }
        else
        {
            inner_index = i - (NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES);
            face_index = inner_index/SCALAR_POINTS_PER_INNER_FACE;
            on_face_index = inner_index - face_index*SCALAR_POINTS_PER_INNER_FACE;
            edge_1 = face_edges[face_index][1];
            edge_2 = face_edges[face_index][2];
            check = 1;
            index_check = -1;
            coord_2 = -1;
            j = 0;
            while (check == 1)
            {
                ++coord_2;
                coord_1_points_amount = POINTS_PER_EDGE - 1 - j;
                index_check = index_check + coord_1_points_amount;
                if (on_face_index <= index_check && on_face_index > index_check - coord_1_points_amount)
                {
                    coord_1 = on_face_index - (index_check - coord_1_points_amount + 1);
                    check = 0;
                }
                ++j;
            }
            reverse_1 = 1;
            reverse_2 = 1;
            if (face_vertices[face_index][1] == edge_vertices[edge_1][0])
                reverse_1 = 0;
            if (face_vertices[face_index][2] == edge_vertices[edge_2][0])
                reverse_2 = 0;
            if (reverse_2 == 0)
                index_1 = NUMBER_OF_PENTAGONS + (edge_2 + 1)*POINTS_PER_EDGE - 1 - coord_2;
            else
                index_1 = NUMBER_OF_PENTAGONS + edge_2*POINTS_PER_EDGE + coord_2;
            if (reverse_1 == 0)
                index_2 = NUMBER_OF_PENTAGONS + edge_1*POINTS_PER_EDGE + coord_2;
            else
                index_2 = NUMBER_OF_PENTAGONS + (edge_1 + 1)*POINTS_PER_EDGE - 1 - coord_2;
            rel_on_line = (1.0 + coord_1)/(1.0 + coord_1_points_amount);
            find_geodetic(latitude_scalar[index_1], longitude_scalar[index_1], latitude_scalar[index_2], longitude_scalar[index_2], rel_on_line, &lat_res, &lon_res);
            latitude_scalar[i] = lat_res;
            longitude_scalar[i] = lon_res;
        }
    }
    free(latitude_ico);
    free(longitude_ico);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
        h_index = i - layer_index*NUMBER_OF_SCALARS_H;
        z_scalar[i] = 0.5*SCALE_HEIGHT*(log(NUMBER_OF_LEVELS/(layer_index + 1.0)) + log(NUMBER_OF_LEVELS/(layer_index + 2.0)));
    }
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
            face_index = (i - NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
            edge_0 = face_edges[face_index][0];
            edge_1 = face_edges[face_index][1];
            edge_2 = face_edges[face_index][2];
            on_face_index = i - (NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
            triangle_on_face_index = on_face_index/3;
            reverse_0 = 1;
            reverse_1 = 1;
            reverse_2 = 1;
            if (face_vertices[face_index][0] == edge_vertices[edge_0][0])
                reverse_0 = 0;
            if (face_vertices[face_index][1] == edge_vertices[edge_1][0])
                reverse_1 = 0;
            if (face_vertices[face_index][2] == edge_vertices[edge_2][0])
                reverse_2 = 0;
            check = 1;
            index_check = -1;
            coord_2 = -1;
            while (check == 1)
            {
                ++coord_2;
                coord_1_points_amount = POINTS_PER_EDGE - coord_2;
                index_check = index_check + coord_1_points_amount;
                if (triangle_on_face_index <= index_check && triangle_on_face_index > index_check - coord_1_points_amount)
                {
                    coord_1 = triangle_on_face_index - (index_check - coord_1_points_amount + 1);
                    check = 0;
                }
            }
            small_triangle_edge_index = on_face_index - 3*triangle_on_face_index;
            if (coord_2 == 0)
                {
                    if (reverse_0 == 0)
                        point_0 = NUMBER_OF_PENTAGONS + edge_0*POINTS_PER_EDGE + coord_1;
                    else
                        point_0 = NUMBER_OF_PENTAGONS + (edge_0 + 1)*POINTS_PER_EDGE - 1 - coord_1;
                }
            else
            {
                point_0 = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES + face_index*SCALAR_POINTS_PER_INNER_FACE + triangle_on_face_index - POINTS_PER_EDGE;
            }
            if (coord_1 == POINTS_PER_EDGE - 1 - coord_2)
                {
                    if (reverse_1 == 0)
                        point_1 = NUMBER_OF_PENTAGONS + edge_1*POINTS_PER_EDGE + coord_2;
                    else
                        point_1 = NUMBER_OF_PENTAGONS + (edge_1 + 1)*POINTS_PER_EDGE - 1 - coord_2;
                }
            else
            {
                point_1 = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES + face_index*SCALAR_POINTS_PER_INNER_FACE + triangle_on_face_index - coord_2;
            }
            if (coord_1 == 0)
                {
                    if (reverse_2 == 0)
                        point_2 = NUMBER_OF_PENTAGONS + (edge_2 + 1)*POINTS_PER_EDGE - 1 - coord_2;
                    else
                        point_2 = NUMBER_OF_PENTAGONS + edge_2*POINTS_PER_EDGE + coord_2;
                }
            else
            {
                point_2 = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES + face_index*SCALAR_POINTS_PER_INNER_FACE + triangle_on_face_index - 1 - coord_2;
            }
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
            dual_scalar_index = face_index*TRIANGLES_PER_FACE + 2*triangle_on_face_index + 1 + coord_2;
            find_center(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_2], longitude_scalar[point_2], &lat_res, &lon_res);
            latitude_scalar_dual[dual_scalar_index] = lat_res;
            longitude_scalar_dual[dual_scalar_index] = lon_res;
            if (coord_2 == 0)
            {
                if (coord_1 == 0)
                    point_3 = face_vertices[face_index][0];
                else
                {
                    if (reverse_0 == 0)
                        point_3 = point_0 - 1;
                    else
                        point_3 = point_0 + 1;
                }
            }
            else if (coord_1 == 0)
            {
                if (reverse_2 == 0)
                    point_3 = point_2 + 1;
                else
                    point_3 = point_2 - 1;
            }
            else
                point_3 = point_0 - 1;
            find_center(latitude_scalar[point_3], longitude_scalar[point_3], latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_2], longitude_scalar[point_2], &lat_res, &lon_res);
            latitude_scalar_dual[dual_scalar_index - 1] = lat_res;
            longitude_scalar_dual[dual_scalar_index - 1] = lon_res;
            if (coord_1 == POINTS_PER_EDGE - 1 - coord_2)
            {
                if (coord_2 == 0)
                    point_3 = face_vertices[face_index][1];
                else
                {
                    if (reverse_1 == 0)
                        point_3 = point_1 - 1;
                    else
                        point_3 = point_1 + 1;
                }
                find_center(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_3], longitude_scalar[point_3], latitude_scalar[point_1], longitude_scalar[point_1], &lat_res, &lon_res);
                latitude_scalar_dual[dual_scalar_index + 1] = lat_res;
                longitude_scalar_dual[dual_scalar_index + 1] = lon_res;
                if (coord_2 == POINTS_PER_EDGE - 1)
                    {
                        point_3 = face_vertices[face_index][2];
                        find_center(latitude_scalar[point_2], longitude_scalar[point_2], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_3], longitude_scalar[point_3], &lat_res, &lon_res);
                        latitude_scalar_dual[dual_scalar_index + 2] = lat_res;
                        longitude_scalar_dual[dual_scalar_index + 2] = lon_res;
                    }
            }
        }
        normal_distance[i] = calculate_distance_h(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], SEMIMAJOR + z_vector[i]);
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        floor_index = layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - floor_index;
        if (h_index >= NUMBER_OF_VECTORS_V)
        {
            gravity[i] = 0;
            h_index -= NUMBER_OF_VECTORS_V;
            if (h_index < NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))
            {
                edge_index = h_index/(POINTS_PER_EDGE + 1);
                on_edge_index = h_index - edge_index*(POINTS_PER_EDGE + 1);
            }
            else
            {
                face_index = (h_index - NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
                on_face_index = h_index - (NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
                triangle_on_face_index = on_face_index/3;
                check = 1;
                index_check = -1;
                coord_2 = -1;
                while (check == 1)
                {
                    ++coord_2;
                    coord_1_points_amount = POINTS_PER_EDGE - coord_2;
                    index_check = index_check + coord_1_points_amount;
                    if (triangle_on_face_index <= index_check && triangle_on_face_index > index_check - coord_1_points_amount)
                    {
                        coord_1 = triangle_on_face_index - (index_check - coord_1_points_amount + 1);
                        check = 0;
                    }
                }
                dual_scalar_index = layer_index*NUMBER_OF_DUAL_SCALARS_H + face_index*TRIANGLES_PER_FACE + 2*triangle_on_face_index + 1 + coord_2;
                z_scalar_dual[dual_scalar_index] = SCALE_HEIGHT*log(NUMBER_OF_LEVELS/(layer_index + 1.0));
                z_scalar_dual[dual_scalar_index - 1] = SCALE_HEIGHT*log(NUMBER_OF_LEVELS/(layer_index + 1.0));
                if (layer_index == NUMBER_OF_LAYERS - 1)
                {
                    z_scalar_dual[dual_scalar_index + NUMBER_OF_DUAL_SCALARS_H] = SCALE_HEIGHT*log(NUMBER_OF_LEVELS/(layer_index + 2.0));
                    z_scalar_dual[dual_scalar_index - 1 + NUMBER_OF_DUAL_SCALARS_H] = SCALE_HEIGHT*log(NUMBER_OF_LEVELS/(layer_index + 2.0));
                }
                if (coord_1 == coord_1_points_amount - 1)
                {
                    z_scalar_dual[dual_scalar_index + 1] = SCALE_HEIGHT*log(NUMBER_OF_LEVELS/(layer_index + 1.0));
                    if (layer_index == NUMBER_OF_LAYERS - 1)
                         z_scalar_dual[dual_scalar_index + 1 + NUMBER_OF_DUAL_SCALARS_H] = SCALE_HEIGHT*log(NUMBER_OF_LEVELS/(layer_index + 2.0));
                    if (coord_2 == POINTS_PER_EDGE - 1)
                    {
                        z_scalar_dual[dual_scalar_index + 2] = SCALE_HEIGHT*log(NUMBER_OF_LEVELS/(layer_index + 1.0));
                        if (layer_index == NUMBER_OF_LAYERS - 1)
                            z_scalar_dual[dual_scalar_index + 2 + NUMBER_OF_DUAL_SCALARS_H] = SCALE_HEIGHT*log(NUMBER_OF_LEVELS/(layer_index + 2.0));
                    }
                }
            }
            z_vector[i] = z_scalar[layer_index*NUMBER_OF_SCALARS_H];
            normal_distance[i] = calculate_distance_h(latitude_scalar[from_index[h_index]], longitude_scalar[from_index[h_index]], latitude_scalar[to_index[h_index]], longitude_scalar[to_index[h_index]], SEMIMAJOR + z_vector[i]);
        }
        else
        {
            gravity[i] = -9.8;
            upper_index = h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H;
            lower_index = h_index + layer_index*NUMBER_OF_SCALARS_H;
            if (layer_index == 0)
                normal_distance[i] = 0;
            else if (layer_index == NUMBER_OF_LAYERS)
                normal_distance[i] = 0;
            else
                normal_distance[i] = z_scalar[upper_index] - z_scalar[lower_index];
            z_vector[i] = SCALE_HEIGHT*log(NUMBER_OF_LEVELS/(1.0 + layer_index));
            if (h_index < NUMBER_OF_PENTAGONS)
                area[i] = 5*4*M_PI*pow(SEMIMAJOR + z_vector[i], 2)/(3*NUMBER_OF_TRIANGLES);
            else
                area[i] = 6*4*M_PI*pow(SEMIMAJOR + z_vector[i], 2)/(3*NUMBER_OF_TRIANGLES);
        }
    }
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
        h_index = i - layer_index*NUMBER_OF_SCALARS_H;
        base_area = area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        radius_1 = SEMIMAJOR + z_vector[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        radius_2 = SEMIMAJOR + z_vector[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER];
        volume[i] = find_volume(base_area, radius_1, radius_2);
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_PER_LAYER; ++i)
    {
        if (i >= NUMBER_OF_DUAL_VECTORS_H)
        {
            latitude_vector_dual[i] = latitude_scalar_dual[i - NUMBER_OF_DUAL_VECTORS_H];
            f_vec[i] = 2*OMEGA*sin(latitude_vector_dual[i]);
        }
        else
        {
            to_index_dual[i] = 0;
            from_index_dual[i] = 1;
            find_geodetic(latitude_scalar_dual[from_index_dual[i]], longitude_scalar_dual[from_index_dual[i]], latitude_scalar_dual[to_index_dual[i]], longitude_scalar_dual[to_index_dual[i]], 0.5, &lat_res, &lon_res);
            latitude_vector_dual[i] = lat_res;
            f_vec[i] = 2*OMEGA*cos(latitude_vector_dual[i])*sin(find_geodetic_direction(latitude_scalar_dual[from_index_dual[i]], longitude_scalar_dual[from_index_dual[i]], latitude_scalar_dual[to_index_dual[i]], longitude_scalar_dual[to_index_dual[i]], 0.5));
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
            z_vector_dual[i] = 0.5*SCALE_HEIGHT*(log(NUMBER_OF_LEVELS/(layer_index + 1.0)) + log(NUMBER_OF_LEVELS/(layer_index + 2.0)));
            normal_distance_dual[i] = z_scalar_dual[upper_index] - z_scalar_dual[lower_index];
            area_dual[i] = 4*M_PI*pow(SEMIMAJOR + z_vector_dual[i], 2)/NUMBER_OF_TRIANGLES;
        }
        else
        {
            z_vector_dual[i] = SCALE_HEIGHT*log(NUMBER_OF_LEVELS/(layer_index + 1.0));
            if (layer_index == 0)
                r_2 = SEMIMAJOR + z_vector_dual[i];
            else
                r_2 = SEMIMAJOR + z_scalar[(layer_index - 1)*NUMBER_OF_SCALARS_H];
            if (layer_index == NUMBER_OF_LAYERS)
                r_1 = SEMIMAJOR;
            else
                r_1 = SEMIMAJOR + z_scalar[layer_index*NUMBER_OF_SCALARS_H];
            primal_vector_index = (NUMBER_OF_LAYERS - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + h_index;
            parallel_distance = normal_distance[primal_vector_index]*(SEMIMAJOR + z_vector_dual[i])/(SEMIMAJOR + z_vector[primal_vector_index]);
            base_distance = parallel_distance*r_1/(SEMIMAJOR + z_vector_dual[i]);
            area_dual[i] = calculate_vertical_face(base_distance, r_1, r_2);
            normal_distance_dual[i] = calculate_distance_h(latitude_scalar_dual[from_index_dual[h_index]], longitude_scalar_dual[from_index_dual[h_index]], latitude_scalar_dual[to_index_dual[h_index]], longitude_scalar_dual[to_index_dual[h_index]], SEMIMAJOR + z_vector_dual[i]);
        }
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_VECTORS_V)
        {
            dual_vector_index = NUMBER_OF_LAYERS*NUMBER_OF_DUAL_VECTORS_PER_LAYER + h_index - NUMBER_OF_VECTORS_V;
            parallel_distance = normal_distance_dual[dual_vector_index]*(SEMIMAJOR + z_vector[i])/(SEMIMAJOR + z_vector_dual[dual_vector_index]);
            r_2 = SEMIMAJOR + z_vector_dual[h_index - NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER];
            r_1 = SEMIMAJOR + z_vector_dual[h_index - NUMBER_OF_VECTORS_V + (layer_index + 1)*NUMBER_OF_DUAL_VECTORS_PER_LAYER];
            base_distance = parallel_distance*r_1/(SEMIMAJOR + z_vector[i]);
            area[i] = calculate_vertical_face(base_distance, r_1, r_2);
        }
    }
    for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
    {
        counter = 0;
        for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
        {
            if (from_index[j] == i || to_index[j] == i)
            {
                adjacent_vector_indices_h[6*i + counter] = j;
                if (from_index[j] == i)
                    adjacent_signs_h[6*i + counter] = 1;
                if (to_index[j] == i)
                    adjacent_signs_h[6*i + counter] = -1;
                counter++;
            }
        }
        if (h_index < NUMBER_OF_PENTAGONS)
        {
            adjacent_vector_indices_h[6*i + 5] = -1;
            adjacent_signs_h[6*i + 5] = 0;
            for (int j = 0; j < 6; ++j)
            {
                if (j < 5)
                {
                    vorticity_indices[6*i + j] = 0;
                    vorticity_signs[6*i + j] = 1;
                }
                else
                {
                    vorticity_indices[6*i + j] = -1;
                    vorticity_signs[6*i + j] = -1;
                }
            }
        }
        else if (h_index < NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES)
        {
            edge_index = (h_index - NUMBER_OF_PENTAGONS)/POINTS_PER_EDGE;
            vertex_index_1 = edge_vertices[edge_index][0];
            vertex_index_2 = edge_vertices[edge_index][1];
            for (int j = 0; j < 6; ++j)
            {
                vorticity_indices[6*i + j] = 0;
                vorticity_signs[6*i + j] = 1;
            }
        }
        else
        {
            inner_index = i - (NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES);
            face_index = inner_index/SCALAR_POINTS_PER_INNER_FACE;
            on_face_index = inner_index - face_index*SCALAR_POINTS_PER_INNER_FACE;
            for (int j = 0; j < 6; ++j)
            {
                vorticity_indices[6*i + j] = 0;
                vorticity_signs[6*i + j] = 1;
            }
        }
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        vector_product_sign[i] = 1;
        for (int j = 0; j < 2; ++j)
        {
            recov_hor_ver_dual_index[2*i + j] = 0;
            recov_hor_ver_dual_weight[2*i + j] = 0;
            recov_hor_par_pri_index[2*i + j] = 0;
            recov_hor_par_pri_weight[2*i + j] = 0;
        }
        for (int j = 0; j < 4; ++j)
        {
            h_curl_indices[4*i + j] = 1;
            h_curl_signs[4*i + j] = 1;
            recov_hor_ver_pri_index[4*i + j] = 0;
            recov_hor_ver_pri_weight[4*i + j] = 0;
        }
        for (int j = 0; j < 11; ++j)
        {
            recov_hor_par_dual_index[11*i + j] = 0;
            recov_hor_par_dual_weight[11*i + j] = 0;
        }
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            recov_ver_1_pri_index[6*i + j] = 0;
            recov_ver_1_pri_weight[6*i + j] = 0;
            recov_ver_1_dual_index[6*i + j] = 0;
            recov_ver_1_dual_weight[6*i + j] = 0;
            recov_ver_2_pri_index[6*i + j] = 0;
            recov_ver_2_pri_weight[6*i + j] = 0;
            recov_ver_2_dual_index[6*i + j] = 0;
            recov_ver_2_dual_weight[6*i + j] = 0;
        }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_H; ++i)
    {
        for (int j = 0; j < 4; ++j)
            {
                h_curl_indices_dual[4*i + j] = 0;
                h_curl_signs_dual[4*i + j] = 1;
            }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_V; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            vorticity_indices_dual[3*i + j] = 0;
            vorticity_signs_dual[3*i + j] = 1;
        }
    }
    int retval, ncid_g_prop;
    if ((retval = nc_create(FILE_NAME, NC_CLOBBER, &ncid_g_prop)))
        ERR(retval);
    int latitude_scalar_id, longitude_scalar_id, latitude_scalar_dual_id, longitude_scalar_dual_id, z_scalar_id, z_vector_id, normal_distance_id, gravity_id, volume_id, area_id, recov_hor_par_dual_weight_id, recov_hor_ver_dual_weight_id, recov_hor_par_pri_weight_id, recov_hor_ver_pri_weight_id, recov_ver_1_pri_weight_id, recov_ver_1_dual_weight_id, recov_ver_2_pri_weight_id, recov_ver_2_dual_weight_id, z_vector_dual_id, normal_distance_dual_id, area_dual_id, f_vec_id, to_index_id, from_index_id, adjacent_vector_indices_h_id, vorticity_indices_id, h_curl_indices_id, recov_hor_par_dual_index_id, recov_hor_ver_dual_index_id, recov_hor_par_pri_index_id, recov_hor_ver_pri_index_id, recov_ver_1_pri_index_id, recov_ver_1_dual_index_id, recov_ver_2_pri_index_id, recov_ver_2_dual_index_id, to_index_dual_id, from_index_dual_id, vorticity_indices_dual_id, h_curl_indices_dual_id, adjacent_signs_h_id, vorticity_signs_id, h_curl_signs_id, vector_product_sign_id, vorticity_signs_dual_id, h_curl_signs_dual_id, vector_index_one_layer_dual;
    int scalar_dimid, scalar_h_dimid, scalar_dual_h_dimid, vector_dimid, scalar_dimid_6, vector_dimid_4, vector_h_dimid, vector_h_dimid_11, vector_h_dimid_2, vector_h_dimid_4, vector_v_dimid_6, vector_dual_dimid, vector_dual_h_dimid, vector_dual_v_dimid_3, vector_dual_h_dimid_4;
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_index", NUMBER_OF_SCALARS, &scalar_dimid)))
        ERR(retval);  
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_h_index", NUMBER_OF_SCALARS_H, &scalar_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_dual_h_index", NUMBER_OF_DUAL_SCALARS_H, &scalar_dual_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_index", NUMBER_OF_VECTORS_H, &vector_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index", NUMBER_OF_VECTORS, &vector_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_6_index", 6*NUMBER_OF_SCALARS_H, &scalar_dimid_6)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_4_index", 4*NUMBER_OF_VECTORS_H, &vector_dimid_4)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_11_index", 11*NUMBER_OF_VECTORS_H, &vector_h_dimid_11)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_2_index", 2*NUMBER_OF_VECTORS_H, &vector_h_dimid_2)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_4_index", 4*NUMBER_OF_VECTORS_H, &vector_h_dimid_4)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_v_6_index", 6*NUMBER_OF_VECTORS_V, &vector_v_dimid_6)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_h_dual", NUMBER_OF_DUAL_VECTORS_H, &vector_dual_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_one_layer_dual", NUMBER_OF_DUAL_VECTORS_PER_LAYER, &vector_index_one_layer_dual)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_dual", NUMBER_OF_DUAL_VECTORS, &vector_dual_dimid)))
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
    if ((retval = nc_def_var(ncid_g_prop, "z_scalar", NC_DOUBLE, 1, &scalar_h_dimid, &z_scalar_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_scalar_id, "units", strlen("m"), "m")))
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
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_par_dual_weight", NC_DOUBLE, 1, &vector_h_dimid_11, &recov_hor_par_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_ver_dual_weight", NC_DOUBLE, 1, &vector_h_dimid_2, &recov_hor_ver_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_par_pri_weight", NC_DOUBLE, 1, &vector_h_dimid_2, &recov_hor_par_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_ver_pri_weight", NC_DOUBLE, 1, &vector_h_dimid_4, &recov_hor_ver_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_1_dual_weight", NC_DOUBLE, 1, &vector_v_dimid_6, &recov_ver_1_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_1_pri_weight", NC_DOUBLE, 1, &vector_v_dimid_6, &recov_ver_1_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_2_pri_weight", NC_DOUBLE, 1, &vector_v_dimid_6, &recov_ver_2_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_2_dual_weight", NC_DOUBLE, 1, &vector_v_dimid_6, &recov_ver_2_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_vector_dual", NC_DOUBLE, 1, &vector_dual_dimid, &z_vector_dual_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_vector_dual_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "normal_distance_dual", NC_DOUBLE, 1, &vector_dual_dimid, &normal_distance_dual_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, normal_distance_dual_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "area_dual", NC_DOUBLE, 1, &vector_dual_dimid, &area_dual_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, area_dual_id, "units", strlen("m^2"), "m^2")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "f_vec", NC_DOUBLE, 1, &vector_index_one_layer_dual, &f_vec_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, f_vec_id, "units", strlen("1/s"), "1/s")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "to_index", NC_LONG, 1, &vector_h_dimid, &to_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "from_index", NC_LONG, 1, &vector_h_dimid, &from_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_vector_indices_h", NC_LONG, 1, &scalar_dimid_6, &adjacent_vector_indices_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_indices", NC_LONG, 1, &scalar_dimid_6, &vorticity_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "h_curl_indices", NC_LONG, 1, &vector_dual_h_dimid_4, &h_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_par_dual_index", NC_LONG, 1, &vector_h_dimid_11, &recov_hor_par_dual_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_ver_dual_index", NC_LONG, 1, &vector_h_dimid_2, &recov_hor_ver_dual_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_par_pri_index", NC_LONG, 1, &vector_h_dimid_2, &recov_hor_par_pri_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_ver_pri_index", NC_LONG, 1, &vector_h_dimid_4, &recov_hor_ver_pri_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_1_pri_index", NC_LONG, 1, &vector_v_dimid_6, &recov_ver_1_pri_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_1_dual_index", NC_LONG, 1, &vector_v_dimid_6, &recov_ver_1_dual_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_2_pri_index", NC_LONG, 1, &vector_v_dimid_6, &recov_ver_2_pri_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_2_dual_index", NC_LONG, 1, &vector_v_dimid_6, &recov_ver_2_dual_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "to_index_dual", NC_LONG, 1, &vector_dual_h_dimid, &to_index_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "from_index_dual", NC_LONG, 1, &vector_dual_h_dimid, &from_index_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_indices_dual", NC_LONG, 1, &vector_dual_v_dimid_3, &vorticity_indices_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "h_curl_indices_dual", NC_LONG, 1, &vector_dual_h_dimid_4, &h_curl_indices_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_signs_h", NC_SHORT, 1, &scalar_dimid_6, &adjacent_signs_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_signs", NC_SHORT, 1, &scalar_dimid_6, &vorticity_signs_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "h_curl_signs", NC_SHORT, 1, &vector_dimid_4, &h_curl_signs_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vector_product_sign", NC_SHORT, 1, &vector_h_dimid, &vector_product_sign_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_signs_dual", NC_SHORT, 1, &vector_dual_v_dimid_3, &vorticity_signs_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "h_curl_signs_dual", NC_SHORT, 1, &vector_dual_h_dimid_4, &h_curl_signs_dual_id)))
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
    if ((retval = nc_put_var_double(ncid_g_prop, recov_hor_par_dual_weight_id, &recov_hor_par_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_hor_ver_dual_weight_id, &recov_hor_ver_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_hor_par_pri_weight_id, &recov_hor_par_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_hor_ver_pri_weight_id, &recov_hor_ver_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_1_pri_weight_id, &recov_ver_1_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_1_dual_weight_id, &recov_ver_1_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_2_pri_weight_id, &recov_ver_2_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_2_dual_weight_id, &recov_ver_2_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, z_vector_dual_id, &z_vector_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, normal_distance_dual_id, &normal_distance_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, area_dual_id, &area_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, f_vec_id, &f_vec[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, to_index_id, &to_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, from_index_id, &from_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_vector_indices_h_id, &adjacent_vector_indices_h[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, vorticity_indices_id, &vorticity_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, h_curl_indices_id, &h_curl_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, recov_hor_par_dual_index_id, &recov_hor_par_dual_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, recov_hor_ver_dual_index_id, &recov_hor_ver_dual_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, recov_hor_par_pri_index_id, &recov_hor_par_pri_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, recov_hor_ver_pri_index_id, &recov_hor_ver_pri_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, recov_ver_1_pri_index_id, &recov_ver_1_pri_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, recov_ver_1_dual_index_id, &recov_ver_1_dual_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, recov_ver_2_pri_index_id, &recov_ver_2_pri_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, recov_ver_2_dual_index_id, &recov_ver_2_dual_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, to_index_dual_id, &to_index_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, from_index_dual_id, &from_index_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, vorticity_indices_dual_id, &vorticity_indices_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, h_curl_indices_dual_id, &h_curl_indices_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_short(ncid_g_prop, adjacent_signs_h_id, &adjacent_signs_h[0])))
        ERR(retval);
    if ((retval = nc_put_var_short(ncid_g_prop, vorticity_signs_id, &vorticity_signs[0])))
        ERR(retval);
    if ((retval = nc_put_var_short(ncid_g_prop, h_curl_signs_id, &h_curl_signs[0])))
        ERR(retval);
    if ((retval = nc_put_var_short(ncid_g_prop, vector_product_sign_id, &vector_product_sign[0])))
        ERR(retval);
    if ((retval = nc_put_var_short(ncid_g_prop, vorticity_signs_dual_id, &vorticity_signs_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_short(ncid_g_prop, h_curl_signs_dual_id, &h_curl_signs_dual[0])))
        ERR(retval);
    if ((retval = nc_close(ncid_g_prop)))
        ERR(retval);
    free(latitude_scalar);
    free(longitude_scalar);
    free(z_scalar);
    free(z_vector);
    free(normal_distance);
    free(gravity);
    free(volume);
    free(area);
    free(recov_hor_par_dual_weight);
    free(recov_hor_ver_dual_weight);
    free(recov_hor_par_pri_weight);
    free(recov_hor_ver_pri_weight);
    free(recov_ver_1_pri_weight);
    free(recov_ver_1_dual_weight);
    free(recov_ver_2_pri_weight);
    free(recov_ver_2_dual_weight);
    free(latitude_scalar_dual);
    free(longitude_scalar_dual);
    free(z_scalar_dual);
    free(latitude_vector_dual);
    free(z_vector_dual);
    free(normal_distance_dual);
    free(f_vec);
    free(to_index);
    free(from_index);
    free(adjacent_vector_indices_h);
    free(vorticity_indices);
    free(h_curl_indices);
    free(recov_hor_par_dual_index);
    free(recov_hor_ver_dual_index);
    free(recov_hor_par_pri_index);
    free(recov_hor_ver_pri_index);
    free(recov_ver_1_pri_index);
    free(recov_ver_1_dual_index);
    free(recov_ver_2_pri_index);
    free(recov_ver_2_dual_index);
    free(to_index_dual);
    free(from_index_dual);
    free(vorticity_indices_dual);
    free(h_curl_indices_dual);
    free(adjacent_signs_h);
    free(vorticity_signs);
    free(h_curl_signs);
    free(vector_product_sign);
    free(vorticity_signs_dual);
    free(h_curl_signs_dual);
    free(area_dual);
    return 0;
}


