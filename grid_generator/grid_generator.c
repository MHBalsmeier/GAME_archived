#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include "/home/max/my_code/game/core/src/enum_and_typedefs.h"
#include "/home/max/my_code/c_source_libs/geos/geos.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define FILE_NAME "grib_files/sphere_res_2_oro_0_geo_prop.nc"

int main(int argc, char *argv[])
{
    const double SCALE_HEIGHT = 8e3;
    const double ATMOS_HEIGHT = SCALE_HEIGHT*log(1 + NUMBER_OF_LAYERS);
    double latitude_ico[12];
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
    double longitude_ico[12];
    longitude_ico[0] = 0;
    longitude_ico[1] = 0;
    longitude_ico[2] = 1*2*M_PI/5;
    longitude_ico[3] = 2*2*M_PI/5;
    longitude_ico[4] = 3*2*M_PI/5;
    longitude_ico[5] = 4*2*M_PI/5;
    longitude_ico[6] = 0;
    longitude_ico[7] = 1*2*M_PI/5 + 2*M_PI/10;
    longitude_ico[8] = 2*2*M_PI/5 + 2*M_PI/10;
    longitude_ico[9] = 3*2*M_PI/5 + 2*M_PI/10;
    longitude_ico[10] = 4*2*M_PI/5 + 2*M_PI/10;
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
    edge_vertices[11][0] = 6;
    edge_vertices[11][1] = 2;
    edge_vertices[12][0] = 2;
    edge_vertices[12][1] = 7;
    edge_vertices[13][0] = 7;
    edge_vertices[13][1] = 3;
    edge_vertices[14][0] = 3;
    edge_vertices[14][1] = 8;
    edge_vertices[15][0] = 8;
    edge_vertices[15][1] = 4;
    edge_vertices[16][0] = 4;
    edge_vertices[16][1] = 9;
    edge_vertices[17][0] = 9;
    edge_vertices[17][1] = 5;
    edge_vertices[18][0] = 5;
    edge_vertices[18][1] = 10;
    edge_vertices[19][0] = 10;
    edge_vertices[19][1] = 1;
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
    int face_vertex_0, face_vertex_1, face_vertex_2, edge_vertex_0, edge_vertex_1;
    for (int i = 0; i < NUMBER_OF_BASIC_TRIANGLES; i++)
    {
        face_vertex_0 = face_vertices[i][0];
        face_vertex_1 = face_vertices[i][1];
        face_vertex_2 = face_vertices[i][2];
        for (int j = 0; j < NUMBER_OF_EDGES; j++)
        {
            edge_vertex_0 = edge_vertices[j][0];
            edge_vertex_1 = edge_vertices[j][1];
            if (edge_vertex_0 == face_vertex_0 && edge_vertex_1 == face_vertex_1 || edge_vertex_0 == face_vertex_1 && edge_vertex_1 == face_vertex_0)
                face_edges[i][0] = j;
            if (edge_vertex_0 == face_vertex_1 && edge_vertex_1 == face_vertex_2 || edge_vertex_0 == face_vertex_2 && edge_vertex_1 == face_vertex_1)
                face_edges[i][1] = j;
            if (edge_vertex_0 == face_vertex_2 && edge_vertex_1 == face_vertex_0 || edge_vertex_0 == face_vertex_0 && edge_vertex_1 == face_vertex_2)
                face_edges[i][2] = j;
        }
    }
    int counter;
    int vertex_edges[12][5];
    for (int i = 0; i < NUMBER_OF_PENTAGONS; i++)
    {
        counter = 0;
        for (int j = 0; j < NUMBER_OF_EDGES; j++)
        {
            edge_vertex_0 = edge_vertices[j][0];
            edge_vertex_1 = edge_vertices[j][1];
                if (edge_vertex_0 == i || edge_vertex_1 == i)
                {
                    vertex_edges[i][counter] = j;
                    counter++;
                }
        }
    }
    int vertex_faces[12][5];
    for (int i = 0; i < NUMBER_OF_PENTAGONS; i++)
    {
        counter = 0;
        for (int j = 0; j < NUMBER_OF_BASIC_TRIANGLES; j++)
        {
            face_vertex_0 = face_vertices[j][0];
            face_vertex_1 = face_vertices[j][1];
            face_vertex_2 = face_vertices[j][2];
                if (face_vertex_0 == i || face_vertex_1 == i || face_vertex_2 == i)
                {
                    vertex_faces[i][counter] = j;
                    counter++;
                }
        }
    }
    double *latitude_scalar, *longitude_scalar, *z_scalar, *latitude_vector, *longitude_vector, *z_vector, *parallel_distance, *normal_distance, *gravity, *volume, *area, *direction, *recov_hor_par_dual_weight, *recov_hor_ver_dual_weight, *recov_hor_par_pri_weight, *recov_hor_ver_pri_weight, *recov_ver_1_pri_weight, *recov_ver_1_dual_weight, *recov_ver_2_pri_weight, *recov_ver_2_dual_weight, *latitude_scalar_dual, *longitude_scalar_dual, *z_scalar_dual, *latitude_vector_dual, *longitude_vector_dual, *z_vector_dual, *normal_distance_dual, *area_dual, *f_vec, *parallel_distance_dual;
    long *to_indices, *from_indices, *adjacent_vector_indices_h, *adjacent_vector_index_lower, *adjacent_vector_index_upper, *adjacent_scalar_indices_h, *adjacent_scalar_index_lower, *adjacent_scalar_index_upper, *vorticity_indices, *h_curl_indices, *recov_hor_par_dual_index, *recov_hor_ver_dual_index, *recov_hor_par_pri_index, *recov_hor_ver_pri_index, *recov_ver_1_pri_index, *recov_ver_1_dual_index, *recov_ver_2_pri_index, *recov_ver_2_dual_index, *to_indices_dual, *from_indices_dual, *adjacent_vector_indices_h_dual, *adjacent_vector_index_upper_dual, *adjacent_vector_index_lower_dual, *adjacent_scalar_indices_h_dual, *adjacent_scalar_index_lower_dual, *adjacent_scalar_index_upper_dual, *vorticity_indices_dual, *h_curl_indices_dual;
    short *adjacent_signs_h, *vorticity_signs, *h_curl_signs, *vector_product_sign, *adjacent_signs_h_dual, *vorticity_signs_dual, *h_curl_signs_dual;
    latitude_scalar = (double *) calloc(NUMBER_OF_SCALARS, sizeof(double));
    longitude_scalar = (double *) calloc(NUMBER_OF_SCALARS, sizeof(double));
    z_scalar = (double *) calloc(NUMBER_OF_SCALARS, sizeof(double));
    latitude_vector = (double *) calloc(NUMBER_OF_VECTORS, sizeof(double));
    longitude_vector = (double *) calloc(NUMBER_OF_VECTORS, sizeof(double));
    z_vector = (double *) calloc(NUMBER_OF_VECTORS, sizeof(double));
    parallel_distance = (double *) calloc(2*NUMBER_OF_VECTORS, sizeof(double));
    normal_distance = (double *) calloc(NUMBER_OF_VECTORS, sizeof(double));
    gravity = (double *) calloc(NUMBER_OF_V_VECTORS, sizeof(double));
    volume = (double *) calloc(NUMBER_OF_SCALARS, sizeof(double));
    area = (double *) calloc(NUMBER_OF_VECTORS, sizeof(double));
    direction = (double *) calloc(NUMBER_OF_VECTORS_H, sizeof(double));
    recov_hor_par_dual_weight = (double *) calloc(11*NUMBER_OF_H_VECTORS, sizeof(double));
    recov_hor_ver_dual_weight = (double *) calloc(2*NUMBER_OF_H_VECTORS, sizeof(double));
    recov_hor_par_pri_weight = (double *) calloc(2*NUMBER_OF_H_VECTORS, sizeof(double));
    recov_hor_ver_pri_weight = (double *) calloc(4*NUMBER_OF_H_VECTORS, sizeof(double));
    recov_ver_1_pri_weight = (double *) calloc(6*NUMBER_OF_V_VECTORS, sizeof(double));
    recov_ver_1_dual_weight = (double *) calloc(6*NUMBER_OF_V_VECTORS, sizeof(double));
    recov_ver_2_pri_weight = (double *) calloc(6*NUMBER_OF_V_VECTORS, sizeof(double));
    recov_ver_2_dual_weight = (double *) calloc(6*NUMBER_OF_V_VECTORS, sizeof(double));
    latitude_scalar_dual = (double *) calloc(NUMBER_OF_DUAL_SCALARS, sizeof(double));
    longitude_scalar_dual = (double *) calloc(NUMBER_OF_DUAL_SCALARS, sizeof(double));
    z_scalar_dual = (double *) calloc(NUMBER_OF_DUAL_SCALARS, sizeof(double));
    latitude_vector_dual = (double *) calloc(NUMBER_OF_DUAL_VECTORS, sizeof(double));
    longitude_vector_dual = (double *) calloc(NUMBER_OF_DUAL_VECTORS, sizeof(double));
    z_vector_dual = (double *) calloc(NUMBER_OF_DUAL_VECTORS, sizeof(double));
    normal_distance_dual = (double *) calloc(NUMBER_OF_DUAL_VECTORS, sizeof(double));
    area_dual = (double *) calloc(NUMBER_OF_DUAL_VECTORS, sizeof(double));
    f_vec = (double *) calloc(NUMBER_OF_DUAL_VECTORS, sizeof(double));
    parallel_distance_dual = (double *) calloc(2*NUMBER_OF_DUAL_VECTORS, sizeof(double));
    recov_ver_2_dual_index = (long *) calloc(6*NUMBER_OF_V_VECTORS, sizeof(long));
    recov_ver_1_dual_index = (long *) calloc(6*NUMBER_OF_V_VECTORS, sizeof(long));
    recov_ver_2_pri_index = (long *) calloc(6*NUMBER_OF_V_VECTORS, sizeof(long));
    recov_hor_ver_pri_index = (long *) calloc(4*NUMBER_OF_H_VECTORS, sizeof(long));
    recov_ver_1_pri_index = (long *) calloc(6*NUMBER_OF_V_VECTORS, sizeof(long));
    recov_hor_ver_dual_index = (long *) calloc(2*NUMBER_OF_H_VECTORS, sizeof(long));
    recov_hor_par_pri_index = (long *) calloc(2*NUMBER_OF_H_VECTORS, sizeof(long));
    to_indices = (long *) calloc(NUMBER_OF_VECTORS, sizeof(long));
    from_indices = (long *) calloc(NUMBER_OF_VECTORS, sizeof(long));
    adjacent_vector_indices_h = (long *) calloc(6*NUMBER_OF_SCALARS, sizeof(long));
    adjacent_vector_index_lower = (long *) calloc(NUMBER_OF_SCALARS, sizeof(long));
    adjacent_vector_index_upper = (long *) calloc(NUMBER_OF_SCALARS, sizeof(long));
    adjacent_scalar_indices_h = (long *) calloc(6*NUMBER_OF_SCALARS, sizeof(long));
    adjacent_scalar_index_lower = (long *) calloc(NUMBER_OF_SCALARS, sizeof(long));
    adjacent_scalar_index_upper = (long *) calloc(NUMBER_OF_SCALARS, sizeof(long));
    vorticity_indices = (long *) calloc(6*NUMBER_OF_SCALARS, sizeof(long));
    h_curl_indices = (long *) calloc(4*NUMBER_OF_VECTORS, sizeof(long));
    recov_hor_par_dual_index = (long *) calloc(11*NUMBER_OF_H_VECTORS, sizeof(long));
    to_indices_dual = (long *) calloc(NUMBER_OF_DUAL_VECTORS, sizeof(long));
    from_indices_dual = (long *) calloc(NUMBER_OF_DUAL_VECTORS, sizeof(long));
    adjacent_vector_indices_h_dual = (long *) calloc(3*NUMBER_OF_DUAL_SCALARS, sizeof(long));
    adjacent_vector_index_upper_dual = (long *) calloc(NUMBER_OF_DUAL_SCALARS, sizeof(long));
    adjacent_vector_index_lower_dual = (long *) calloc(NUMBER_OF_DUAL_SCALARS, sizeof(long));
    adjacent_scalar_indices_h_dual = (long *) calloc(3*NUMBER_OF_DUAL_SCALARS, sizeof(long));
    adjacent_scalar_index_lower_dual = (long *) calloc(NUMBER_OF_DUAL_SCALARS, sizeof(long));
    adjacent_scalar_index_upper_dual = (long *) calloc(NUMBER_OF_DUAL_SCALARS, sizeof(long));
    vorticity_indices_dual = (long *) calloc(3*NUMBER_OF_DUAL_V_VECTORS, sizeof(long));
    h_curl_indices_dual = (long *) calloc(4*NUMBER_OF_DUAL_H_VECTORS, sizeof(long));
    adjacent_signs_h = (short *) calloc(6*NUMBER_OF_SCALARS, sizeof(short));
    vorticity_signs = (short *) calloc(6*NUMBER_OF_SCALARS, sizeof(short));
    h_curl_signs = (short *) calloc(4*NUMBER_OF_VECTORS, sizeof(short));
    vector_product_sign = (short *) calloc(NUMBER_OF_H_VECTORS, sizeof(short));
    adjacent_signs_h_dual = (short *) calloc(3*NUMBER_OF_DUAL_SCALARS, sizeof(short));
    vorticity_signs_dual = (short *) calloc(3*NUMBER_OF_DUAL_V_VECTORS, sizeof(short));
    h_curl_signs_dual = (short *) calloc(4*NUMBER_OF_DUAL_H_VECTORS, sizeof(short));
    double sigma, sigma_1, sigma_2, lat_edge_1, lon_edge_1, lat_edge_2, lon_edge_2, point_frac, lat_1, lon_1, lat_2, lon_2, rel_on_line, lat_res, lon_res, base_area, radius_1, radius_2, z_1, z_2;
    long h_index, edge_index, face_index, on_face_index, layer_index, level_index, inner_index, check, edge_0, edge_1, edge_2, edge_index_1, edge_index_2, coord_1_points_amount, index_check, j, coord_1, coord_2, vertex_index_0, vertex_index_1, vertex_index_2, reverse_0, reverse_1, reverse_2, index_1, index_2, index_3, points_right, small_triangle_edge_index, triangle_on_face_index, vert_index, floor_index, on_edge_index, point_0, point_1, point_2, dual_scalar_index, vector_index;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
        h_index = i - layer_index*NUMBER_OF_SCALARS_H;
        sigma = 0.5*(SCALE_HEIGHT/ATMOS_HEIGHT)*(log((1.0 + NUMBER_OF_LAYERS)/(layer_index + 1.0)) + log((1.0 + NUMBER_OF_LAYERS)/(layer_index + 2.0)));
        z_scalar[i] = ATMOS_HEIGHT*sigma;
        if (layer_index == 0)
            adjacent_vector_index_upper[i] = -1;
        else
            adjacent_vector_index_upper[i] = h_index + layer_index*(NUMBER_OF_VECTORS_H + NUMBER_OF_SCALARS_H);
        if (layer_index == NUMBER_OF_LAYERS - 1)
            adjacent_vector_index_lower[i] = -1;
        else
            adjacent_vector_index_lower[i] = h_index + (layer_index + 1)*(NUMBER_OF_VECTORS_H + NUMBER_OF_SCALARS_H);
        if (h_index < NUMBER_OF_PENTAGONS)
        {
            latitude_scalar[i] = latitude_ico[h_index];
            longitude_scalar[i] = longitude_ico[h_index];
        }
        else if (h_index < NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES)
        {
            edge_index = (h_index - NUMBER_OF_PENTAGONS)/POINTS_PER_EDGE;
            point_frac = (1.0 + h_index - (NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE))/(POINTS_PER_EDGE + 1.0);
            vertex_index_0 = edge_vertices[edge_index][0];
            vertex_index_1 = edge_vertices[edge_index][1];
            find_geodetic(latitude_ico[vertex_index_0], longitude_ico[vertex_index_0], latitude_ico[vertex_index_1], longitude_ico[vertex_index_1], point_frac, &lat_res, &lon_res);
            latitude_scalar[i] = lat_res;
            longitude_scalar[i] = lon_res;
        }
        else
        {
            inner_index = i - layer_index*NUMBER_OF_SCALARS_H - (NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES);
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
                index_1 = (edge_2 + 1)*POINTS_PER_EDGE - 1 - coord_2 + NUMBER_OF_PENTAGONS;
            else
                index_1 = edge_2*POINTS_PER_EDGE  + coord_2 + NUMBER_OF_PENTAGONS;
            if (reverse_1 == 0)
                index_2 = edge_1*POINTS_PER_EDGE + coord_2 + NUMBER_OF_PENTAGONS;
            else
                index_2 = (edge_1 + 1)*POINTS_PER_EDGE - 1 - coord_1 + NUMBER_OF_PENTAGONS;
            rel_on_line = (1.0 + coord_1)/(coord_1_points_amount + 1.0);
            find_geodetic(latitude_scalar[index_1], longitude_scalar[index_1], latitude_scalar[index_2], longitude_scalar[index_2], rel_on_line, &lat_res, &lon_res);
            latitude_scalar[i] = lat_res;
            longitude_scalar[i] = lon_res;
        }
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        vert_index = i/(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H);
        floor_index = vert_index*(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H);
        h_index = i - floor_index;
        if (h_index >= NUMBER_OF_SCALARS_H)
        {
            h_index -= NUMBER_OF_SCALARS_H;
            if (h_index < NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))
            {
                edge_index = h_index/POINTS_PER_EDGE;
                on_edge_index = h_index - edge_index*POINTS_PER_EDGE;
                if(on_edge_index == 0)
                {
                    from_indices[i] = edge_vertices[edge_index][0];
                    to_indices[i] = NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE;
                }
                else if (on_edge_index == POINTS_PER_EDGE - 1)
                {
                    from_indices[i] = NUMBER_OF_PENTAGONS + (edge_index + 1)*POINTS_PER_EDGE - 1;
                    to_indices[i] = edge_vertices[edge_index][1];
                }
                else
                {
                    from_indices[i] = NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE + on_edge_index - 1;
                    to_indices[i] = NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE + on_edge_index;
                }
            }
            else
            {
                face_index = (h_index - NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
                edge_0 = face_edges[face_index][0];
                edge_1 = face_edges[face_index][1];
                edge_2 = face_edges[face_index][2];
                on_face_index = h_index - (NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
                triangle_on_face_index = on_face_index/3;
                reverse_0, reverse_1, reverse_2 = 1;
                if (face_vertices[face_index][0] == edge_vertices[edge_0][0])
                    reverse_1 = 0;
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
                        coord_1 = on_face_index - (index_check - coord_1_points_amount + 1);
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
                    point_0 = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES + face_index*SCALAR_POINTS_PER_INNER_FACE + triangle_on_face_index - POINTS_PER_EDGE - coord_2 + 1;
                if (coord_1 == POINTS_PER_EDGE - 1 - coord_2)
                    {
                        if (reverse_1 == 0)
                            point_1 = NUMBER_OF_PENTAGONS + edge_1*POINTS_PER_EDGE + coord_2;
                        else
                            point_1 = NUMBER_OF_PENTAGONS + (edge_1 + 1)*POINTS_PER_EDGE - 1 - coord_1;
                    }
                else
                    point_1 = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES + face_index*SCALAR_POINTS_PER_INNER_FACE + triangle_on_face_index - coord_2;
                if (coord_1 == 0)
                    {
                        if (reverse_2 == 0)
                            point_2 = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*(edge_2 + 1) - 1 - coord_2;
                        else
                            point_2 = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*edge_2 + coord_2;
                    }
                else
                    point_2 = NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES + face_index*SCALAR_POINTS_PER_INNER_FACE + triangle_on_face_index - 1 - coord_2;
                if (small_triangle_edge_index == 0)
                {
                    from_indices[i] = point_0;
                    to_indices[i] = point_2;
                }
                if (small_triangle_edge_index == 1)
                {
                    from_indices[i] = point_0;
                    to_indices[i] = point_1;
                }
                if (small_triangle_edge_index == 2)
                {
                    from_indices[i] = point_2;
                    to_indices[i] = point_1;
                }
                dual_scalar_index = triangle_on_face_index + face_index*TRIANGLES_PER_FACE + vert_index*NUMBER_OF_DUAL_SCALARS_H;
                find_center(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_2], longitude_scalar[point_2], &lat_res, &lon_res);
                latitude_scalar_dual[dual_scalar_index] = lat_res;
                longitude_scalar_dual[dual_scalar_index] = lon_res;
                if (vert_index == 0)
                    adjacent_vector_index_upper_dual[dual_scalar_index] = -1;
                else
                    adjacent_vector_index_upper_dual[dual_scalar_index] = dual_scalar_index + (vert_index - 1)*NUMBER_OF_DUAL_VECTORS_H;
                if (vert_index == NUMBER_OF_LAYERS)
                    adjacent_vector_index_lower_dual[dual_scalar_index] = -1;
                else
                    adjacent_vector_index_lower_dual[dual_scalar_index] = dual_scalar_index + vert_index*NUMBER_OF_DUAL_VECTORS_H;
                z_scalar_dual[dual_scalar_index] = SCALE_HEIGHT/ATMOS_HEIGHT*log(NUMBER_OF_LEVELS/(vert_index + 1));
            }
            find_geodetic(latitude_scalar[from_indices[i]], longitude_scalar[from_indices[i]], latitude_scalar[to_indices[i]], longitude_scalar[to_indices[i]], 0.5, &lat_res, &lon_res);
            latitude_vector[i] = lat_res;
            longitude_vector[i] = lon_res;
            z_vector[i] = z_scalar[vert_index*NUMBER_OF_SCALARS_H];
            if (layer_index == 0)
                direction[h_index] = find_geodetic_direction(latitude_scalar[from_indices[i]], longitude_scalar[from_indices[i]], latitude_scalar[to_indices[i]], longitude_scalar[to_indices[i]], 0.5);
            normal_distance[i] = calculate_distance_h(latitude_vector[from_indices[i]], longitude_vector[from_indices[i]], latitude_vector[to_indices[i]], longitude_vector[to_indices[i]], SEMIMAJOR + z_vector[i]);
            parallel_distance[2*i + 0] = 1;
            parallel_distance[2*i + 1] = 1;
            area[i] = calculate_vertical_face(parallel_distance[2*i + 0], SEMIMAJOR + ATMOS_HEIGHT*sigma_2, SEMIMAJOR + ATMOS_HEIGHT*sigma_1);
        }
        else
        {
            if (vert_index == 0)
            {
                from_indices[i] = h_index + (vert_index + 1)*NUMBER_OF_SCALARS_H;
                to_indices[i] = -1;
                latitude_vector[i] = latitude_scalar[from_indices[i]];
                longitude_vector[i] = longitude_scalar[from_indices[i]];
                normal_distance[i] = -1;
            }
            else if (vert_index == NUMBER_OF_LAYERS)
            {
                from_indices[i] = -1;
                to_indices[i] = h_index + vert_index*NUMBER_OF_SCALARS_H;
                latitude_vector[i] = latitude_scalar[to_indices[i]];
                longitude_vector[i] = longitude_scalar[to_indices[i]];
                normal_distance[i] = -1;
            }
            else
            {
                from_indices[i] = h_index + vert_index*NUMBER_OF_SCALARS_H;
                to_indices[i] = h_index + (vert_index - 1)*NUMBER_OF_SCALARS_H;
                latitude_vector[i] = latitude_scalar[to_indices[i]];
                longitude_vector[i] = longitude_scalar[to_indices[i]];
                normal_distance[i] = z_scalar[to_indices[i]] - z_scalar[from_indices[i]];
            }
            z_vector[i] = ATMOS_HEIGHT*(SCALE_HEIGHT/ATMOS_HEIGHT)*log((1.0 + NUMBER_OF_LAYERS)/(1.0 + vert_index));
            parallel_distance[2*i + 0] = 1;
            parallel_distance[2*i + 1] = 1;
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
        base_area = area[adjacent_vector_index_lower[i]];
        radius_2 = SEMIMAJOR + z_vector[adjacent_vector_index_upper[i]];
        radius_1 = SEMIMAJOR + z_vector[adjacent_vector_index_lower[i]];
        volume[i] = find_volume(base_area, radius_1, radius_2);
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_DUAL_SCALARS_H);
        h_index = i - layer_index*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_DUAL_SCALARS_H);
        if (h_index >= NUMBER_OF_DUAL_VECTORS_H)
        {
            z_vector_dual[i] = ATMOS_HEIGHT*0.5*(SCALE_HEIGHT/ATMOS_HEIGHT)*(log(NUMBER_OF_LEVELS/(layer_index + 1.0)) + log(NUMBER_OF_LEVELS/(layer_index + 2.0)));
            to_indices_dual[i] = h_index - NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_SCALARS_H;
            from_indices_dual[i] = h_index  - NUMBER_OF_DUAL_VECTORS_H + (layer_index + 1)*NUMBER_OF_DUAL_SCALARS_H;
            latitude_vector_dual[i] = latitude_scalar_dual[to_indices_dual[i]];
            longitude_vector_dual[i] = longitude_scalar_dual[to_indices_dual[i]];
            normal_distance_dual[i] = z_scalar_dual[to_indices_dual[i]] - z_scalar_dual[from_indices_dual[i]];
            area_dual[i] = 4*M_PI*pow(SEMIMAJOR + z_vector_dual[i], 2)/NUMBER_OF_TRIANGLES;
            parallel_distance_dual[2*i + 0] = 1;
            parallel_distance_dual[2*i + 1] = 1;
        }
        else
        {
            z_vector_dual[i] = ATMOS_HEIGHT*(SCALE_HEIGHT/ATMOS_HEIGHT)*log(NUMBER_OF_LEVELS/(layer_index + 1.0));
            parallel_distance_dual[2*i + 0] = 1;
            parallel_distance_dual[2*i + 1] = 1;
            area_dual[i] = 0;
            to_indices_dual[i] = 0;
            from_indices_dual[i] = 0;
            normal_distance_dual[i] = 0;
            latitude_vector_dual[i] = 0;
            longitude_vector_dual[i] = 0;
        }
        if (h_index >= NUMBER_OF_DUAL_VECTORS_H)
            f_vec[i] = 2*OMEGA*sin(latitude_vector_dual[i]);
        else
            f_vec[i] = 2*OMEGA*cos(latitude_vector_dual[i])*sin(find_geodetic_direction(latitude_scalar_dual[from_indices_dual[i]], longitude_scalar_dual[from_indices_dual[i]], latitude_scalar_dual[to_indices_dual[i]], longitude_scalar_dual[to_indices_dual[i]], 0.5));
    }
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
        if (layer_index == 0)
            adjacent_scalar_index_upper[i] = -1;
        else
            adjacent_scalar_index_upper[i] = i - NUMBER_OF_SCALARS_H;
        if (layer_index == NUMBER_OF_LAYERS - 1)
            adjacent_scalar_index_lower[i] = -1;
        else
            adjacent_scalar_index_lower[i] = i + NUMBER_OF_SCALARS_H;
        counter = 0;
        for (int j = 0; j < NUMBER_OF_VECTORS_H; j++)
        {
            vector_index = layer_index*(NUMBER_OF_VECTORS_H + NUMBER_OF_SCALARS_H) + NUMBER_OF_SCALARS_H + j;
            if (from_indices[vector_index] == i || to_indices[vector_index] == i)
            {
                adjacent_vector_indices_h[6*i + counter] = vector_index;
                if (from_indices[vector_index] == i)
                {
                    adjacent_signs_h[6*i + counter] = 1;
                    adjacent_scalar_indices_h[6*i + counter] = to_indices[vector_index];
                }
                if (to_indices[vector_index] == i)
                {
                    adjacent_signs_h[6*i + counter] = -1;
                    adjacent_scalar_indices_h[6*i + counter] = from_indices[vector_index];
                }
                counter++;
            }
        }
        if (h_index < NUMBER_OF_PENTAGONS)
        {
            adjacent_scalar_indices_h[6*i + 5] = -1;
            adjacent_vector_indices_h[6*i + 5] = -1;
            adjacent_signs_h[6*i + 5] = 0;
            for (int j = 0; j < 6; j++)
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
            for (int j = 0; j < 6; j++)
            {
                vorticity_indices[6*i + j] = 0;
                vorticity_signs[6*i + j] = 1;
            }
        }
        else
        {
            inner_index = i - layer_index*NUMBER_OF_SCALARS_H - (NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES);
            face_index = inner_index/SCALAR_POINTS_PER_INNER_FACE;
            on_face_index = inner_index - face_index*SCALAR_POINTS_PER_INNER_FACE;
            for (int j = 0; j < 6; j++)
            {
                vorticity_indices[6*i + j] = 0;
                vorticity_signs[6*i + j] = 1;
            }
        }
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        for (int j = 0; j < 4; j++)
            {
                h_curl_indices[4*i + j] = 1;
                h_curl_signs[4*i + j] = 1;
            }
        }
    for (int i = 0; i < NUMBER_OF_H_VECTORS; ++i)
        {
        vector_product_sign[i] = 1;
        for (int j = 0; j < 2; j++)
        {
            recov_hor_ver_dual_index[2*i + j] = 0;
            recov_hor_ver_dual_weight[2*i + j] = 0;
            recov_hor_par_pri_index[2*i + j] = 0;
            recov_hor_par_pri_weight[2*i + j] = 0;
        }
        for (int j = 0; j < 4; j++)
        {
            h_curl_indices[4*i + j] = 1;
            h_curl_signs[4*i + j] = 1;
            recov_hor_ver_pri_index[4*i + j] = 0;
            recov_hor_ver_pri_weight[4*i + j] = 0;
        }
        for (int j = 0; j < 11; j++)
        {
            recov_hor_par_dual_index[11*i + j] = 0;
            recov_hor_par_dual_weight[11*i + j] = 0;
        }
    }
    for (int i = 0; i < NUMBER_OF_V_VECTORS; ++i)
    {
        gravity[i] = -9.8;
        for (int j = 0; j < 6; j++)
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
    for (int i = 0; i < NUMBER_OF_DUAL_SCALARS; ++i)
    {
        level_index = i/NUMBER_OF_DUAL_SCALARS_H;
        h_index = i - level_index*NUMBER_OF_DUAL_SCALARS_H;
        face_index = h_index/TRIANGLES_PER_FACE;
        on_face_index =  h_index - face_index*TRIANGLES_PER_FACE;
        edge_1 = face_edges[face_index][0];
        edge_1 = face_edges[face_index][1];
        edge_2 = face_edges[face_index][2];
        if (level_index == 0)
            adjacent_scalar_index_upper_dual[i] = -1;
        else
            adjacent_scalar_index_upper_dual[i] = i - NUMBER_OF_DUAL_SCALARS_H;
        if (level_index == NUMBER_OF_LEVELS - 1)
            adjacent_scalar_index_lower_dual[i] = -1;
        else
            adjacent_scalar_index_lower_dual[i] = i + NUMBER_OF_DUAL_SCALARS_H;
        counter = 0;
        for (int j = 0; j < NUMBER_OF_DUAL_VECTORS_H; j++)
        {
            vector_index = level_index*(NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_DUAL_SCALARS_H) + j;
            if (from_indices_dual[vector_index] == i || to_indices_dual[vector_index] == i)
            {
                adjacent_vector_indices_h_dual[3*i + counter] = vector_index;
                if (from_indices_dual[vector_index] == i)
                {
                    adjacent_signs_h_dual[3*i + counter] = 1;
                    adjacent_scalar_indices_h_dual[3*i + counter] = to_indices_dual[vector_index];
                }
                if (to_indices_dual[vector_index] == i)
                {
                    adjacent_signs_h_dual[3*i + counter] = -1;
                    adjacent_scalar_indices_h_dual[3*i + counter] = from_indices_dual[vector_index];
                }
                counter++;
            }
        }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_H_VECTORS; ++i)
    {
        for (int j = 0; j < 4; j++)
            {
                h_curl_indices_dual[4*i + j] = 0;
                h_curl_signs_dual[4*i + j] = 1;
            }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_V_VECTORS; ++i)
    {
        for (int j = 0; j < 3; j++)
        {
            vorticity_indices_dual[3*i + j] = 0;
            vorticity_signs_dual[3*i + j] = 1;
        }
    }
    int retval, ncid_g_prop;
    if ((retval = nc_create(FILE_NAME, NC_CLOBBER, &ncid_g_prop)))
        ERR(retval);
    int latitude_scalar_id, longitude_scalar_id, z_scalar_id, latitude_vector_id, longitude_vector_id, z_vector_id, parallel_distance_id, normal_distance_id, gravity_id, volume_id, area_id, recov_hor_par_dual_weight_id, recov_hor_ver_dual_weight_id, recov_hor_par_pri_weight_id, recov_hor_ver_pri_weight_id, recov_ver_1_pri_weight_id, recov_ver_1_dual_weight_id, recov_ver_2_pri_weight_id, recov_ver_2_dual_weight_id, latitude_scalar_dual_id, longitude_scalar_dual_id, z_scalar_dual_id, latitude_vector_dual_id, longitude_vector_dual_id, z_vector_dual_id, normal_distance_dual_id, area_dual_id, f_vec_id, parallel_distance_dual_id, to_indices_id, from_indices_id, adjacent_vector_indices_h_id, adjacent_vector_index_lower_id, adjacent_vector_index_upper_id, adjacent_scalar_indices_h_id, adjacent_scalar_index_lower_id, adjacent_scalar_index_upper_id, vorticity_indices_id, h_curl_indices_id, recov_hor_par_dual_index_id, recov_hor_ver_dual_index_id, recov_hor_par_pri_index_id, recov_hor_ver_pri_index_id, recov_ver_1_pri_index_id, recov_ver_1_dual_index_id, recov_ver_2_pri_index_id, recov_ver_2_dual_index_id, to_indices_dual_id, from_indices_dual_id, adjacent_vector_indices_h_dual_id, adjacent_vector_index_upper_dual_id, adjacent_vector_index_lower_dual_id, adjacent_scalar_indices_h_dual_id, adjacent_scalar_index_lower_dual_id, adjacent_scalar_index_upper_dual_id, vorticity_indices_dual_id, h_curl_indices_dual_id, adjacent_signs_h_id, vorticity_signs_id, h_curl_signs_id, vector_product_sign_id, adjacent_signs_h_dual_id, vorticity_signs_dual_id, h_curl_signs_dual_id;
    int scalar_dimid, vector_dimid, vector_v_dimid, vector_dimid_2, scalar_dimid_6, vector_dimid_4, vector_h_dimid, vector_h_dimid_11, vector_h_dimid_2, vector_h_dimid_4, vector_v_dimid_6, scalar_dual_dimid, vector_dual_dimid, vector_dual_h_dimid, vector_dual_dimid_2, scalar_dual_dimid_3, vector_dual_v_dimid_3, vector_dual_h_dimid_4;
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_index", NUMBER_OF_SCALARS, &scalar_dimid)))
        ERR(retval);  
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index", NUMBER_OF_VECTORS, &vector_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_2_index", 2*NUMBER_OF_VECTORS, &vector_dimid_2)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_6_index", 6*NUMBER_OF_SCALARS, &scalar_dimid_6)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_4_index", 4*NUMBER_OF_VECTORS, &vector_dimid_4)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_index", NUMBER_OF_VECTORS_H, &vector_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_v_index", NUMBER_OF_V_VECTORS, &vector_v_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_11_index", 11*NUMBER_OF_H_VECTORS, &vector_h_dimid_11)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_2_index", 2*NUMBER_OF_H_VECTORS, &vector_h_dimid_2)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_4_index", 4*NUMBER_OF_H_VECTORS, &vector_h_dimid_4)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_v_6_index", 6*NUMBER_OF_V_VECTORS, &vector_v_dimid_6)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_index_dual", NUMBER_OF_DUAL_SCALARS, &scalar_dual_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_h_dual", NUMBER_OF_DUAL_VECTORS_H, &vector_dual_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_dual", NUMBER_OF_DUAL_VECTORS, &vector_dual_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_dual_2_index", 2*NUMBER_OF_DUAL_VECTORS, &vector_dual_dimid_2)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_dual_3_index", 3*NUMBER_OF_DUAL_SCALARS, &scalar_dual_dimid_3)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_dual_v_3_index", 3*NUMBER_OF_DUAL_V_VECTORS, &vector_dual_v_dimid_3)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_dual_h_4_index", 4*NUMBER_OF_DUAL_H_VECTORS, &vector_dual_h_dimid_4)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "latitude_scalar", NC_DOUBLE, 1, &scalar_dimid, &latitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "longitude_scalar", NC_DOUBLE, 1, &scalar_dimid, &longitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_scalar", NC_DOUBLE, 1, &scalar_dimid, &z_scalar_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_scalar_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "latitude_vector", NC_DOUBLE, 1, &vector_dimid, &latitude_vector_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "longitude_vector", NC_DOUBLE, 1, &vector_dimid, &longitude_vector_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_vector", NC_DOUBLE, 1, &vector_dimid, &z_vector_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, z_vector_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "parallel_distance", NC_DOUBLE, 1, &vector_dimid_2, &parallel_distance_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, parallel_distance_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "normal_distance", NC_DOUBLE, 1, &vector_dimid, &normal_distance_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, normal_distance_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "gravity", NC_DOUBLE, 1, &vector_v_dimid, &gravity_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "latitude_scalar_dual", NC_DOUBLE, 1, &scalar_dual_dimid, &latitude_scalar_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "longitude_scalar_dual", NC_DOUBLE, 1, &scalar_dual_dimid, &longitude_scalar_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "z_scalar_dual", NC_DOUBLE, 1, &scalar_dual_dimid, &z_scalar_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "latitude_vector_dual", NC_DOUBLE, 1, &vector_dual_dimid, &latitude_vector_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "longitude_vector_dual", NC_DOUBLE, 1, &vector_dual_dimid, &longitude_vector_dual_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "f_vec", NC_DOUBLE, 1, &vector_dual_dimid, &f_vec_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, f_vec_id, "units", strlen("1/s"), "1/s")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "parallel_distance_dual", NC_DOUBLE, 1, &vector_dual_dimid_2, &parallel_distance_dual_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, parallel_distance_dual_id, "units", strlen("m"), "m")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "to_indices", NC_LONG, 1, &vector_dual_dimid, &to_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "from_indices", NC_LONG, 1, &vector_dual_dimid, &from_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_vector_indices_h", NC_LONG, 1, &scalar_dimid_6, &adjacent_vector_indices_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_vector_index_lower", NC_LONG, 1, &scalar_dimid, &adjacent_vector_index_lower_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_vector_index_upper", NC_LONG, 1, &scalar_dimid, &adjacent_vector_index_upper_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_scalar_indices_h", NC_LONG, 1, &scalar_dimid_6, &adjacent_scalar_indices_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_scalar_index_lower", NC_LONG, 1, &scalar_dimid, &adjacent_scalar_index_lower_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_scalar_index_upper", NC_LONG, 1, &scalar_dimid, &adjacent_scalar_index_upper_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "to_indices_dual", NC_LONG, 1, &vector_dual_dimid, &to_indices_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "from_indices_dual", NC_LONG, 1, &vector_dual_dimid, &from_indices_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_vector_indices_h_dual", NC_LONG, 1, &scalar_dual_dimid_3, &adjacent_vector_indices_h_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_vector_index_upper_dual", NC_LONG, 1, &scalar_dual_dimid, &adjacent_vector_index_upper_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_vector_index_lower_dual", NC_LONG, 1, &scalar_dual_dimid, &adjacent_vector_index_lower_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_scalar_indices_h_dual", NC_LONG, 1, &scalar_dual_dimid_3, &adjacent_scalar_indices_h_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_scalar_index_lower_dual", NC_LONG, 1, &scalar_dual_dimid, &adjacent_scalar_index_lower_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_scalar_index_upper_dual", NC_LONG, 1, &scalar_dual_dimid, &adjacent_scalar_index_upper_dual_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_signs_h_dual", NC_SHORT, 1, &scalar_dual_dimid_3, &adjacent_signs_h_dual_id)))
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
    if ((retval = nc_put_var_double(ncid_g_prop, z_scalar_id, &z_scalar[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, latitude_vector_id, &latitude_vector[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, longitude_vector_id, &longitude_vector[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, z_vector_id, &z_vector[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, parallel_distance_id, &parallel_distance[0])))
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
    if ((retval = nc_put_var_double(ncid_g_prop, latitude_scalar_dual_id, &latitude_scalar_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, longitude_scalar_dual_id, &longitude_scalar_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, z_scalar_dual_id, &z_scalar_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, latitude_vector_dual_id, &latitude_vector_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, longitude_vector_dual_id, &longitude_vector_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, z_vector_dual_id, &z_vector_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, normal_distance_dual_id, &normal_distance_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, area_dual_id, &area_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, f_vec_id, &f_vec[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, parallel_distance_dual_id, &parallel_distance_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, to_indices_id, &to_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, from_indices_id, &from_indices[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_vector_indices_h_id, &adjacent_vector_indices_h[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_vector_index_lower_id, &adjacent_vector_index_lower[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_vector_index_upper_id, &adjacent_vector_index_upper[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_scalar_indices_h_id, &adjacent_scalar_indices_h[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_scalar_index_lower_id, &adjacent_scalar_index_lower[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_scalar_index_upper_id, &adjacent_scalar_index_upper[0])))
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
    if ((retval = nc_put_var_long(ncid_g_prop, to_indices_dual_id, &to_indices_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, from_indices_dual_id, &from_indices_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_vector_indices_h_dual_id, &adjacent_vector_indices_h_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_vector_index_upper_dual_id, &adjacent_vector_index_upper_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_vector_index_lower_dual_id, &adjacent_vector_index_lower_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_scalar_indices_h_dual_id, &adjacent_scalar_indices_h_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_scalar_index_lower_dual_id, &adjacent_scalar_index_lower_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, adjacent_scalar_index_upper_dual_id, &adjacent_scalar_index_upper_dual[0])))
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
    if ((retval = nc_put_var_short(ncid_g_prop, adjacent_signs_h_dual_id, &adjacent_signs_h_dual[0])))
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
    free(latitude_vector);
    free(longitude_vector);
    free(z_vector);
    free(parallel_distance);
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
    free(longitude_vector_dual);
    free(z_vector_dual);
    free(normal_distance_dual);
    free(area_dual);
    free(f_vec);
    free(parallel_distance_dual);
    free(to_indices);
    free(from_indices);
    free(adjacent_vector_indices_h);
    free(adjacent_vector_index_lower);
    free(adjacent_vector_index_upper);
    free(adjacent_scalar_indices_h);
    free(adjacent_scalar_index_lower);
    free(adjacent_scalar_index_upper);
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
    free(to_indices_dual);
    free(from_indices_dual);
    free(adjacent_vector_indices_h_dual);
    free(adjacent_vector_index_upper_dual);
    free(adjacent_vector_index_lower_dual);
    free(adjacent_scalar_indices_h_dual);
    free(adjacent_scalar_index_lower_dual);
    free(adjacent_scalar_index_upper_dual);
    free(vorticity_indices_dual);
    free(h_curl_indices_dual);
    free(adjacent_signs_h);
    free(vorticity_signs);
    free(h_curl_signs);
    free(vector_product_sign);
    free(adjacent_signs_h_dual);
    free(vorticity_signs_dual);
    free(h_curl_signs_dual);
    return 0;
}


