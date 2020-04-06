#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <math.h>
#include "/lib/geos/include/geos.h"
#include "/lib/conv/include/conv.h"
#include "/lib/indextools/include/index_tools.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define OMEGA (7.292115e-5)

const short MODE = 2;
const double TOA = 30000.0;
const double SCALE_HEIGHT = 8000.0;
const short ORO_ID = 0;
const double ORTH_CRITERION_DEG = 89.8;

enum grid_integers {
RES_ID = 4,
NUMBER_OF_BASIC_TRIANGLES = 20,
NUMBER_OF_PENTAGONS = 12,
NUMBER_OF_HEXAGONS = (int) (10*(pow(2, 2*RES_ID) - 1)),
NUMBER_OF_EDGES = 3*NUMBER_OF_BASIC_TRIANGLES/2,
NUMBER_OF_LAYERS = 6,
NUMBER_OF_LEVELS = NUMBER_OF_LAYERS + 1,
NUMBER_OF_SCALARS_H = NUMBER_OF_PENTAGONS + NUMBER_OF_HEXAGONS,
NUMBER_OF_VECTORS_H = (5*NUMBER_OF_PENTAGONS/2 + 6/2*NUMBER_OF_HEXAGONS),
NUMBER_OF_VECTORS_V = NUMBER_OF_SCALARS_H,
NUMBER_OF_H_VECTORS = NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_H,
NUMBER_OF_V_VECTORS = NUMBER_OF_LEVELS*NUMBER_OF_VECTORS_V,
NUMBER_OF_VECTORS_PER_LAYER = NUMBER_OF_VECTORS_H + NUMBER_OF_VECTORS_V,
NUMBER_OF_TRIANGLES = (int) (NUMBER_OF_BASIC_TRIANGLES*(pow(4, RES_ID))),
NUMBER_OF_SCALARS = NUMBER_OF_SCALARS_H*NUMBER_OF_LAYERS,
NUMBER_OF_VECTORS = NUMBER_OF_H_VECTORS + NUMBER_OF_V_VECTORS,
NUMBER_OF_DUAL_SCALARS_H = NUMBER_OF_TRIANGLES,
NUMBER_OF_DUAL_VECTORS_H = 3*NUMBER_OF_TRIANGLES/2,
NUMBER_OF_DUAL_VECTORS_V = NUMBER_OF_DUAL_SCALARS_H,
NUMBER_OF_DUAL_H_VECTORS = NUMBER_OF_LEVELS*NUMBER_OF_DUAL_VECTORS_H,
NUMBER_OF_DUAL_V_VECTORS = NUMBER_OF_LAYERS*NUMBER_OF_DUAL_VECTORS_V,
NUMBER_OF_DUAL_VECTORS_PER_LAYER = NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_DUAL_VECTORS_V,
NUMBER_OF_DUAL_SCALARS = NUMBER_OF_LEVELS*NUMBER_OF_DUAL_SCALARS_H,
NUMBER_OF_DUAL_VECTORS = NUMBER_OF_DUAL_H_VECTORS + NUMBER_OF_DUAL_V_VECTORS,
TRIANGLES_PER_FACE = NUMBER_OF_TRIANGLES/NUMBER_OF_BASIC_TRIANGLES,
POINTS_PER_EDGE = (int) (pow(2, RES_ID) - 1),
SCALAR_POINTS_PER_INNER_FACE = (int) (0.5*(pow(2, RES_ID) - 2)*(pow(2, RES_ID) - 1)),
VECTOR_POINTS_PER_INNER_FACE = (int) (1.5*(pow(2, RES_ID) - 1)*pow(2, RES_ID))};

int find_angle_change(double, double, double *);
int find_coords_from_triangle_on_face_index(long, long, long *, long *, long *);
int find_triangle_on_face_index_from_coords(long, long, long, long *);
int find_sigma_from_level(int, double *);
int find_sigma_from_layer(int, double *);
int find_triangle_indices_from_h_vector_index(int, long, long *, long *, long *, long *, long *, long *, long *, short *, short[][3], short[][3], short[][2], short[][3]);
int find_triangle_edge_points(long, short, long, long *, long *, long *, long *, long *, long *, long *, short[][3], short[][3], short[][3]);
int find_points_per_edge(long, long *);
int find_scalar_points_per_inner_face(long, long *);
int upscale_scalar_point(long, long, long *);
int write_scalar_coordinates(long, long, long, long, long, long, short, double[], double[], double[], double[], double[]);
int find_triangles_per_face(long, long *);
int find_triangle_edge_points_from_dual_scalar_on_face_index(long, short, long, long *, long *, long *, short[][3], short[][3], short[][3]);
int find_triangle_on_face_index_from_dual_scalar_on_face_index(long, long, long *, short *, short *, short *);

int main(int argc, char *argv[])
{
    short OUTPUT_FILE_LENGTH = 100;
    char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE_PRE, "nc_files/B%dL%dT%d_M%d_O%d", RES_ID, NUMBER_OF_LAYERS, (int) TOA, MODE, ORO_ID);
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "nc_files/B%dL%dT%d_M%d_O%d.nc", RES_ID, NUMBER_OF_LAYERS, (int) TOA, MODE, ORO_ID);
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
    short edge_vertices[NUMBER_OF_EDGES][2];
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
    short face_vertices[20][3];
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
    short face_edges[20][3];
    short face_edges_reverse[20][3];
    short edge_other_vertex_index, check_index;
    for (int i = 0; i < 20; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 30; ++k)
            {
                if (edge_vertices[k][0] == face_vertices[i][j] || edge_vertices[k][1] == face_vertices[i][j])
                {
                    if (edge_vertices[k][0] == face_vertices[i][j])
                        edge_other_vertex_index = 1;
                    if (edge_vertices[k][1] == face_vertices[i][j])
                        edge_other_vertex_index = 0;
                    if (j == 0)
                        check_index = 1;
                    if (j == 1)
                        check_index = 2;
                    if (j == 2)
                        check_index = 0;
                    if (edge_vertices[k][edge_other_vertex_index] == face_vertices[i][check_index])
                    {
                        face_edges[i][j] = k;
                        if (edge_other_vertex_index == 1)
                            face_edges_reverse[i][j] = 0;
                        if (edge_other_vertex_index == 0)
                            face_edges_reverse[i][j] = 1;
                    }
                }
            }
        }
    }
    double *x_unity = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *y_unity = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *z_unity = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *latitude_scalar = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *z_scalar = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *z_vector = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *normal_distance = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *latitude_vector = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *longitude_vector = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *direction = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *gravity = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *volume = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *area = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *recov_hor_par_pri_weight = malloc(4*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_hor_ver_pri_weight = malloc(4*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_hor_par_dual_weight = malloc(2*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_hor_ver_dual_weight = malloc(2*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_ver_0_pri_weight = malloc(12*NUMBER_OF_VECTORS_V*sizeof(double));
    double *recov_ver_1_pri_weight = malloc(12*NUMBER_OF_VECTORS_V*sizeof(double));
    double *recov_ver_0_dual_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *recov_ver_1_dual_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *latitude_scalar_dual = malloc(NUMBER_OF_DUAL_SCALARS_H*sizeof(double));
    double *longitude_scalar_dual = malloc(NUMBER_OF_DUAL_SCALARS_H*sizeof(double));
    double *z_scalar_dual = malloc(NUMBER_OF_DUAL_SCALARS*sizeof(double));
    double *latitude_vector_dual = malloc(NUMBER_OF_DUAL_VECTORS_PER_LAYER*sizeof(double));
    double *z_vector_dual = malloc(NUMBER_OF_DUAL_VECTORS*sizeof(double));
    double *normal_distance_dual = malloc(NUMBER_OF_DUAL_VECTORS*sizeof(double));
    double *direction_dual = malloc(NUMBER_OF_DUAL_VECTORS_H*sizeof(double));
    double *area_dual = malloc(NUMBER_OF_DUAL_VECTORS*sizeof(double));
    double *f_vec = malloc(NUMBER_OF_DUAL_VECTORS_PER_LAYER*sizeof(double));
    double *triangle_face_unit_sphere = malloc(NUMBER_OF_DUAL_VECTORS_V*sizeof(double));
    double *pent_hex_face_unity_sphere = malloc(NUMBER_OF_VECTORS_V*sizeof(double));
    long *to_index = malloc(NUMBER_OF_VECTORS_H*sizeof(long));
    long *from_index = malloc(NUMBER_OF_VECTORS_H*sizeof(long));
    long *recov_hor_par_pri_index = malloc(4*NUMBER_OF_VECTORS_H*sizeof(long));
    long *recov_hor_ver_pri_index = malloc(4*NUMBER_OF_VECTORS_H*sizeof(long));
    long *recov_hor_par_dual_index = malloc(2*NUMBER_OF_VECTORS_H*sizeof(long));
    long *recov_hor_ver_dual_index = malloc(2*NUMBER_OF_VECTORS_H*sizeof(long));
    long *recov_ver_0_pri_index = malloc(12*NUMBER_OF_VECTORS_V*sizeof(long));
    long *recov_ver_1_pri_index = malloc(12*NUMBER_OF_VECTORS_V*sizeof(long));
    long *recov_ver_0_dual_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(long));
    long *recov_ver_1_dual_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(long));
    long *adjacent_vector_indices_h = malloc(6*NUMBER_OF_SCALARS_H*sizeof(long));
    long *vorticity_indices = malloc(3*NUMBER_OF_DUAL_VECTORS_V*sizeof(long));
    long *h_curl_indices = malloc(4*NUMBER_OF_DUAL_VECTORS_H*sizeof(long));
    long *to_index_dual = malloc(NUMBER_OF_DUAL_VECTORS_H*sizeof(long));
    long *from_index_dual = malloc(NUMBER_OF_DUAL_VECTORS_H*sizeof(long));
    long *vorticity_indices_dual = malloc(6*NUMBER_OF_VECTORS_V*sizeof(long));
    long *h_curl_indices_dual = malloc(4*NUMBER_OF_VECTORS_H*sizeof(long));
    short *adjacent_signs_h = malloc(6*NUMBER_OF_SCALARS_H*sizeof(short));
    short *vorticity_signs = malloc(3*NUMBER_OF_DUAL_VECTORS_V*sizeof(short));
    short *h_curl_signs = malloc(4*NUMBER_OF_DUAL_VECTORS_H*sizeof(short));
    short *vorticity_signs_dual = malloc(6*NUMBER_OF_VECTORS_V*sizeof(short));
    short *h_curl_signs_dual = malloc(4*NUMBER_OF_VECTORS_H*sizeof(short));
    double lat_edge_1, lon_edge_1, lat_edge_2, lon_edge_2, lat_1, lon_1, lat_2, lon_2, lat_res, lon_res, base_area, base_distance, z_0, z_1, radius_0, radius_1, parallel_distance, direction_change, x_point_0, y_point_0, z_point_0, x_point_1, y_point_1, z_point_1, x_res, y_res, z_res, sigma_z, sigma_z_dual_scalar, rel_on_line;
    int face_index, face_index_0, face_index_1, retval;
    long h_index, on_face_index, layer_index, level_index, inner_index, upper_index, lower_index, coord_0_points_amount, j, coord_0, coord_1, index_0, index_1, index_2, points_right, triangle_on_face_index, on_edge_index, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_index, counter, primal_vector_index, dual_vector_index, number_of_triangles_per_face, edgepoint_0, edgepoint_1, edgepoint_2, edgepoint_3, points_per_edge, dual_scalar_on_face_index, base_index_old, base_index_down_triangles, base_index_up_triangles, old_triangle_on_line_index;
    short first_face_found, edge_rel_to_face_0, edge_rel_to_face_1, edge_index, sign, small_triangle_edge_index, dump, points_upwards, points_downwards, last_triangle_bool;
    if (NUMBER_OF_VECTORS_H != NUMBER_OF_DUAL_VECTORS_H)
        printf("It is NUMBER_OF_VECTORS_H != NUMBER_OF_DUAL_VECTORS_H.\n");
    long test_index;
    for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
    {
        upscale_scalar_point(RES_ID, i, &test_index);
        if (test_index != i)
            printf("problem with upscale_scalar_point detected\n");
    }
    if (MODE == 0)
    {
        for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
        {
            if (i < NUMBER_OF_PENTAGONS)
            {
                latitude_scalar[i] = latitude_ico[i];
                longitude_scalar[i] = longitude_ico[i];
                retval = find_global_normal(latitude_ico[i], longitude_ico[i], &x_res, &y_res, &z_res);
                x_unity[i] = x_res;
                y_unity[i] = y_res;
                z_unity[i] = z_res;
            }
            else if (i < NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES)
            {
                edge_index = (i - NUMBER_OF_PENTAGONS)/POINTS_PER_EDGE;
                rel_on_line = (1.0 + i - (NUMBER_OF_PENTAGONS + edge_index*POINTS_PER_EDGE))/(1.0 + POINTS_PER_EDGE);
                retval = find_between_point(x_unity[edge_vertices[edge_index][0]], y_unity[edge_vertices[edge_index][0]], z_unity[edge_vertices[edge_index][0]], x_unity[edge_vertices[edge_index][1]], y_unity[edge_vertices[edge_index][1]], z_unity[edge_vertices[edge_index][1]], rel_on_line, &x_res, &y_res, &z_res);
                x_unity[i] = x_res;
                y_unity[i] = y_res;
                z_unity[i] = z_res;
                retval = find_geos(x_res, y_res, z_res, &lat_res, &lon_res);
                latitude_scalar[i] = lat_res;
                longitude_scalar[i] = lon_res;
            }
            else
            {
                inner_index = i - (NUMBER_OF_PENTAGONS + POINTS_PER_EDGE*NUMBER_OF_EDGES);
                face_index = inner_index/SCALAR_POINTS_PER_INNER_FACE;
                on_face_index = inner_index - face_index*SCALAR_POINTS_PER_INNER_FACE;
                retval = find_coords_from_triangle_on_face_index(on_face_index + POINTS_PER_EDGE, RES_ID, &coord_0, &coord_1, &coord_0_points_amount);
                if (face_edges_reverse[face_index][2] == 0)
                    index_0 = NUMBER_OF_PENTAGONS + (face_edges[face_index][2] + 1)*POINTS_PER_EDGE - coord_1;
                else
                    index_0 = NUMBER_OF_PENTAGONS + face_edges[face_index][2]*POINTS_PER_EDGE + (coord_1 - 1);
                if (face_edges_reverse[face_index][1] == 0)
                    index_1 = NUMBER_OF_PENTAGONS + face_edges[face_index][1]*POINTS_PER_EDGE + (coord_1 - 1);
                else
                    index_1 = NUMBER_OF_PENTAGONS + (face_edges[face_index][1] + 1)*POINTS_PER_EDGE - coord_1;
                rel_on_line = (1.0 + coord_0)/(1.0 + coord_0_points_amount);
                retval = find_global_normal(latitude_scalar[index_0], longitude_scalar[index_0], &x_point_0, &y_point_0, &z_point_0);
                retval = find_global_normal(latitude_scalar[index_1], longitude_scalar[index_1], &x_point_1, &y_point_1, &z_point_1);
                retval = find_between_point(x_point_0, y_point_0, z_point_0, x_point_1, y_point_1, z_point_1, rel_on_line, &x_res, &y_res, &z_res);
                x_unity[i] = x_res;
                y_unity[i] = y_res;
                z_unity[i] = z_res;
                retval = find_geos(x_res, y_res, z_res, &lat_res, &lon_res);
                latitude_scalar[i] = lat_res;
                longitude_scalar[i] = lon_res;
            }
        }
    }
    if (MODE == 1 || MODE == 2)
    {
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
    }
    free(latitude_ico);
    free(longitude_ico);
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
        retval = find_sigma_from_layer(layer_index, &sigma_z);
        z_scalar[i] = TOA*sigma_z;
        if (z_scalar[i] <= 0)
            printf("z_scalar contains a non-positive value.\n");
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
            if (MODE == 2)
                retval += find_voronoi_center_sphere(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_2], longitude_scalar[point_2], &lat_res, &lon_res);
            else
            {
                retval = find_barycenter_cart(x_unity[point_0], y_unity[point_0], z_unity[point_0], x_unity[point_1], y_unity[point_1], z_unity[point_1], x_unity[point_2], y_unity[point_2], z_unity[point_2], &x_res, &y_res, &z_res);
                retval = find_geos(x_res, y_res, z_res, &lat_res, &lon_res);
            }
            latitude_scalar_dual[dual_scalar_index] = lat_res;
            longitude_scalar_dual[dual_scalar_index] = lon_res;
            retval += calc_triangle_face(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_2], longitude_scalar[point_2], &triangle_face);
            triangle_face_unit_sphere[dual_scalar_index] = triangle_face;
            if (MODE == 2)
                 retval += find_voronoi_center_sphere(latitude_scalar[point_3], longitude_scalar[point_3], latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_2], longitude_scalar[point_2], &lat_res, &lon_res);
            else
            {
                retval = find_barycenter_cart(x_unity[point_3], y_unity[point_3], z_unity[point_3], x_unity[point_0], y_unity[point_0], z_unity[point_0], x_unity[point_2], y_unity[point_2], z_unity[point_2], &x_res, &y_res, &z_res);
                retval = find_geos(x_res, y_res, z_res, &lat_res, &lon_res);
            }
            latitude_scalar_dual[dual_scalar_index - 1] = lat_res;
            longitude_scalar_dual[dual_scalar_index - 1] = lon_res;
            retval += calc_triangle_face(latitude_scalar[point_3], longitude_scalar[point_3], latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_2], longitude_scalar[point_2], &triangle_face);
            triangle_face_unit_sphere[dual_scalar_index - 1] = triangle_face;
            if (coord_0 == coord_0_points_amount - 1)
            {
                if (MODE == 2)
                    retval += find_voronoi_center_sphere(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_4], longitude_scalar[point_4], latitude_scalar[point_1], longitude_scalar[point_1], &lat_res, &lon_res);
                else
                {
                    retval = find_barycenter_cart(x_unity[point_0], y_unity[point_0], z_unity[point_0], x_unity[point_4], y_unity[point_4], z_unity[point_4], x_unity[point_1], y_unity[point_1], z_unity[point_1], &x_res, &y_res, &z_res);
                    retval = find_geos(x_res, y_res, z_res, &lat_res, &lon_res);
                }
                latitude_scalar_dual[dual_scalar_index + 1] = lat_res;
                longitude_scalar_dual[dual_scalar_index + 1] = lon_res;
                retval += calc_triangle_face(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_4], longitude_scalar[point_4], latitude_scalar[point_1], longitude_scalar[point_1], &triangle_face);
                triangle_face_unit_sphere[dual_scalar_index + 1] = triangle_face;
                if (coord_1 == POINTS_PER_EDGE - 1)
                {
                    if (MODE == 2)
                        retval += find_voronoi_center_sphere(latitude_scalar[point_2], longitude_scalar[point_2], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_5], longitude_scalar[point_5], &lat_res, &lon_res);
                    else
                    {
                        retval = find_barycenter_cart(x_unity[point_2], y_unity[point_2], z_unity[point_2], x_unity[point_1], y_unity[point_1], z_unity[point_1], x_unity[point_5], y_unity[point_5], z_unity[point_5], &x_res, &y_res, &z_res);
                        retval = find_geos(x_res, y_res, z_res, &lat_res, &lon_res);
                    }
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
            printf("triangle_face_unit_sphere contains a non-positive value.\n");
        if (fabs(triangle_face_unit_sphere[i]/triangle_avg_unit_sphere_ideal - 1) > 0.4)
            printf("Triangles on unit sphere have significantly different surfaces.\n");
    }
    if (fabs(triangle_sum_unit_sphere/(4*M_PI) - 1) > 1e-13)
        printf("Sum of faces of triangles on unit sphere does not match face of unit sphere.\n");
    free(x_unity);
    free(y_unity);
    free(z_unity);
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_PER_LAYER; ++i)
    {
        if (i >= NUMBER_OF_DUAL_VECTORS_H)
        {
            latitude_vector_dual[i] = latitude_scalar_dual[i - NUMBER_OF_DUAL_VECTORS_H];
            f_vec[i] = 2*OMEGA*sin(latitude_vector_dual[i]);
        }
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
            retval = find_min_dist_rel_on_line(latitude_scalar_dual[from_index_dual[i]], longitude_scalar_dual[from_index_dual[i]], latitude_scalar_dual[to_index_dual[i]], longitude_scalar_dual[to_index_dual[i]], latitude_vector_dual[i], longitude_vector[i], &rel_on_line);
            if (fabs(rel_on_line - 0.5) > 0.14)
                printf("Bisection warning.\n");
            direction_dual[i] = find_geodetic_direction(latitude_scalar_dual[from_index_dual[i]], longitude_scalar_dual[from_index_dual[i]], latitude_scalar_dual[to_index_dual[i]], longitude_scalar_dual[to_index_dual[i]], rel_on_line);
            f_vec[i] = 2*OMEGA*cos(latitude_vector_dual[i])*sin(direction_dual[i]);
        }
    }
    short trouble_detected = 0;
    for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
    {
        counter = 0;
        for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
        {
            if (from_index[j] == i || to_index[j] == i)
            {
                adjacent_vector_indices_h[6*i + counter] = j;
                vorticity_indices_dual[6*i + counter] = j;
                sign = 1;
                if (from_index[j] == i)
                {
                    adjacent_signs_h[6*i + counter] = 1;
                    find_angle_change(direction[j], direction_dual[j], &direction_change);
                    if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
                        sign = -1;
                }
                if (to_index[j] == i)
                {
                    adjacent_signs_h[6*i + counter] = -1;
                    find_angle_change(direction[j], direction_dual[j], &direction_change);
                    if (rad2deg(direction_change) > ORTH_CRITERION_DEG)
                        sign = -1;
                }
                vorticity_signs_dual[6*i + counter] = sign;
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
            printf("Trouble detected, place 1.\n");
        if (i < NUMBER_OF_PENTAGONS)
        {
            adjacent_vector_indices_h[6*i + 5] = -1;
            adjacent_signs_h[6*i + 5] = 0;
            vorticity_indices_dual[6*i + 5] = 0;
            vorticity_signs_dual[6*i + 5] = 0;
        }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_V; ++i)
    {
        counter = 0;
        for (int j = 0; j < NUMBER_OF_DUAL_VECTORS_H; ++j)
        {
            if (from_index_dual[j] == i || to_index_dual[j] == i)
            {
                vorticity_indices[3*i + counter] = j;
                sign = 1;
                if (from_index_dual[j] == i)
                {
                    find_angle_change(direction_dual[j], direction[j], &direction_change);
                    if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
                        sign = -1;
                }
                if (to_index_dual[j] == i)
                {
                    find_angle_change(direction_dual[j], direction[j], &direction_change);
                    if (rad2deg(direction_change) > ORTH_CRITERION_DEG)
                        sign = -1;
                }
                vorticity_signs[3*i + counter] = sign;
                ++counter;
            }
        }
        if (counter != 3)
            printf("Trouble detected, place 0.\n");
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        find_angle_change(direction[i], direction_dual[i], &direction_change);
        if (fabs(rad2deg(direction_change)) < ORTH_CRITERION_DEG || fabs(rad2deg(direction_change)) > 90 + (90 - ORTH_CRITERION_DEG))
            printf("grid non-orthogonal\n");
    }
    short check_0, check_1, check_2;
    for (int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
    {
        if (i < NUMBER_OF_PENTAGONS)
        {
            double *lat_points = malloc(5*sizeof(double));
            double *lon_points = malloc(5*sizeof(double));
            long *cell_vector_indices = malloc(5*sizeof(long));
            for (int j = 0; j < 5; ++j)
                cell_vector_indices[j] = adjacent_vector_indices_h[6*i + j];
            counter = 0;
            for (int j = 0; j < NUMBER_OF_DUAL_VECTORS_V; ++j)
            {
                retval = in_bool_calculator_long(cell_vector_indices, 5, vorticity_indices[3*j + 0], &check_0);
                retval = in_bool_calculator_long(cell_vector_indices, 5, vorticity_indices[3*j + 1], &check_1);
                retval = in_bool_calculator_long(cell_vector_indices, 5, vorticity_indices[3*j + 2], &check_2);
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
            long *cell_vector_indices = malloc(6*sizeof(long));
            for (int j = 0; j < 6; ++j)
                cell_vector_indices[j] = adjacent_vector_indices_h[6*i + j];
            counter = 0;
            for (int j = 0; j < NUMBER_OF_DUAL_VECTORS_V; ++j)
            {
                retval = in_bool_calculator_long(cell_vector_indices, 6, vorticity_indices[3*j + 0], &check_0);
                retval = in_bool_calculator_long(cell_vector_indices, 6, vorticity_indices[3*j + 1], &check_1);
                retval = in_bool_calculator_long(cell_vector_indices, 6, vorticity_indices[3*j + 2], &check_2);
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
            printf("pent_hex_face_unity_sphere contains a non-positive value.\n");
        if (fabs(pent_hex_face_unity_sphere[i]/pent_hex_avg_unity_sphere_ideal - 1) > 0.4)
            printf("Pentagons and hexagons on unity sphere have significantly different surfaces.\n");
    }
    if (fabs(pent_hex_sum_unity_sphere/(4*M_PI) - 1) > 1e-12)
        printf("Sum of faces of pentagons and hexagons on unity sphere does not match face of unit sphere.\n");
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
                dual_scalar_index = layer_index*NUMBER_OF_DUAL_SCALARS_H + face_index*TRIANGLES_PER_FACE + 1 + 2*triangle_on_face_index + coord_1;
                retval = find_sigma_from_level(layer_index, &sigma_z_dual_scalar);
                z_scalar_dual[dual_scalar_index] = TOA*sigma_z_dual_scalar;
                z_scalar_dual[dual_scalar_index - 1] = TOA*sigma_z_dual_scalar;
                if (layer_index == NUMBER_OF_LAYERS - 1)
                {
                    z_scalar_dual[dual_scalar_index + NUMBER_OF_DUAL_SCALARS_H] = 0;
                    z_scalar_dual[dual_scalar_index - 1 + NUMBER_OF_DUAL_SCALARS_H] = 0;
                }
                if (coord_0 == coord_0_points_amount - 1)
                {
                    z_scalar_dual[dual_scalar_index + 1] = TOA*sigma_z_dual_scalar;
                    if (layer_index == NUMBER_OF_LAYERS - 1)
                         z_scalar_dual[dual_scalar_index + 1 + NUMBER_OF_DUAL_SCALARS_H] = 0;
                    if (coord_1 == POINTS_PER_EDGE - 1)
                    {
                        z_scalar_dual[dual_scalar_index + 2] = TOA*sigma_z_dual_scalar;
                        if (layer_index == NUMBER_OF_LAYERS - 1)
                            z_scalar_dual[dual_scalar_index + 2 + NUMBER_OF_DUAL_SCALARS_H] = 0;
                    }
                }
            }
            gravity[i] = 0;
            z_vector[i] = z_scalar[layer_index*NUMBER_OF_SCALARS_H];
            if (z_vector[i] <= 0)
                printf("z_vector contains a non-positive value at a horizontal grid point.\n");
            normal_distance[i] = calculate_distance_h(latitude_scalar[from_index[h_index - NUMBER_OF_VECTORS_V]], longitude_scalar[from_index[h_index - NUMBER_OF_VECTORS_V]], latitude_scalar[to_index[h_index - NUMBER_OF_VECTORS_V]], longitude_scalar[to_index[h_index - NUMBER_OF_VECTORS_V]], SEMIMAJOR + z_vector[i]);
        }
        else
        {
            gravity[i] = -9.80616;
            upper_index = h_index + (layer_index - 1)*NUMBER_OF_SCALARS_H;
            lower_index = h_index + layer_index*NUMBER_OF_SCALARS_H;
            retval = find_sigma_from_level(layer_index, &sigma_z);
            if (layer_index == 0)
                normal_distance[i] = TOA*(1 - (layer_index + 0.0)/NUMBER_OF_LAYERS) - z_scalar[lower_index];
            else if (layer_index == NUMBER_OF_LAYERS)
                normal_distance[i] = z_scalar[upper_index] - 0;
            else
                normal_distance[i] = z_scalar[upper_index] - z_scalar[lower_index];
            z_vector[i] = TOA*sigma_z;
            if (z_vector[i] < 0)
                printf("z_vector contains a negative value.\n");
            area[i] = pent_hex_face_unity_sphere[h_index]*pow(SEMIMAJOR + z_vector[i], 2);
        }
    }
    free(pent_hex_face_unity_sphere);
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        if (normal_distance[i] <= 0)
            printf("normal_distance contains a non-positive value.\n");
    }
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        layer_index = i/NUMBER_OF_SCALARS_H;
        h_index = i - layer_index*NUMBER_OF_SCALARS_H;
        base_area = area[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        radius_0 = SEMIMAJOR + z_vector[h_index + (layer_index + 1)*NUMBER_OF_VECTORS_PER_LAYER];
        radius_1 = SEMIMAJOR + z_vector[h_index + layer_index*NUMBER_OF_VECTORS_PER_LAYER];
        volume[i] = find_volume(base_area, radius_0, radius_1);
    }
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        if (volume[i] <= 0)
            printf("volume contains a non-positive value.\n");
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_DUAL_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_DUAL_VECTORS_H)
        {
            retval = find_sigma_from_layer(layer_index, &sigma_z);
            upper_index = h_index - NUMBER_OF_DUAL_VECTORS_H + layer_index*NUMBER_OF_DUAL_SCALARS_H;
            lower_index = h_index - NUMBER_OF_DUAL_VECTORS_H + (layer_index + 1)*NUMBER_OF_DUAL_SCALARS_H;
            z_vector_dual[i] = TOA*sigma_z;
            normal_distance_dual[i] = z_scalar_dual[upper_index] - z_scalar_dual[lower_index];
            area_dual[i] = pow(SEMIMAJOR + z_vector_dual[i], 2)*triangle_face_unit_sphere[h_index - NUMBER_OF_DUAL_VECTORS_H];
        }
        else
        {
            retval = find_sigma_from_level(layer_index, &sigma_z);
            z_vector_dual[i] = TOA*sigma_z;
            if (layer_index == 0)
                radius_1 = SEMIMAJOR + z_vector_dual[i];
            else
                radius_1 = SEMIMAJOR + z_scalar[(layer_index - 1)*NUMBER_OF_SCALARS_H];
            if (layer_index == NUMBER_OF_LAYERS)
                radius_0 = SEMIMAJOR;
            else
                radius_0 = SEMIMAJOR + z_scalar[layer_index*NUMBER_OF_SCALARS_H];
            primal_vector_index = (NUMBER_OF_LAYERS - 1)*NUMBER_OF_VECTORS_PER_LAYER + NUMBER_OF_VECTORS_V + h_index;
            parallel_distance = normal_distance[primal_vector_index]*(SEMIMAJOR + z_vector_dual[i])/(SEMIMAJOR + z_vector[primal_vector_index]);
            base_distance = parallel_distance*radius_0/(SEMIMAJOR + z_vector_dual[i]);
            area_dual[i] = calculate_vertical_face(base_distance, radius_0, radius_1);
            normal_distance_dual[i] = calculate_distance_h(latitude_scalar_dual[from_index_dual[h_index]], longitude_scalar_dual[from_index_dual[h_index]], latitude_scalar_dual[to_index_dual[h_index]], longitude_scalar_dual[to_index_dual[h_index]], SEMIMAJOR + z_vector_dual[i]);
        }
    }
    free(triangle_face_unit_sphere);
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; ++i)
    {
        if (area_dual[i] <= 0)
            printf("area_dual contains a non-positive value.\n");
        if (normal_distance_dual[i] <= 0)
            printf("normal_distance_dual contains a non-positive value.\n");
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        layer_index = i/NUMBER_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NUMBER_OF_VECTORS_PER_LAYER;
        if (h_index >= NUMBER_OF_VECTORS_V)
        {
            dual_vector_index = NUMBER_OF_LAYERS*NUMBER_OF_DUAL_VECTORS_PER_LAYER + h_index - NUMBER_OF_VECTORS_V;
            parallel_distance = normal_distance_dual[dual_vector_index]*(SEMIMAJOR + z_vector[i])/(SEMIMAJOR + z_vector_dual[dual_vector_index]);
            radius_1 = SEMIMAJOR + z_vector_dual[h_index - NUMBER_OF_VECTORS_V + layer_index*NUMBER_OF_DUAL_VECTORS_PER_LAYER];
            radius_0 = SEMIMAJOR + z_vector_dual[h_index - NUMBER_OF_VECTORS_V + (layer_index + 1)*NUMBER_OF_DUAL_VECTORS_PER_LAYER];
            base_distance = parallel_distance*radius_0/(SEMIMAJOR + z_vector[i]);
            area[i] = calculate_vertical_face(base_distance, radius_0, radius_1);
        }
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        if (area[i] <= 0)
            printf("area contains a non-positive value.\n");
    }
    long *adjacent_scalar_indices_for_cross = malloc(6*sizeof(long));
    short *face_of_cell_indices = malloc(2*sizeof(short));
    short bool_0, bool_1, first_found;
    long cell_0_for_cross, cell_1_for_cross;
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            if (j == 0)
                recov_hor_ver_dual_index[2*i + j] = to_index_dual[i];
            else
                recov_hor_ver_dual_index[2*i + j] = from_index_dual[i];
            recov_hor_ver_dual_weight[2*i + j] = 0.5;
            sign = 1;
            recov_hor_par_dual_index[2*i + j] = i + j*NUMBER_OF_DUAL_VECTORS_PER_LAYER;
            find_angle_change(direction[i], direction_dual[i], &direction_change);
            if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
                sign = -1;
            recov_hor_par_dual_weight[2*i + j] = sign*0.5;
        }
        first_found = 0;
        for (int j = 0; j < NUMBER_OF_SCALARS_H; ++j)
        {
            for (int k = 0; k < 5; ++k)
            {
                if (to_index[adjacent_vector_indices_h[6*j + k]] != j)
                    adjacent_scalar_indices_for_cross[k] = to_index[adjacent_vector_indices_h[6*j + k]];
                else
                    adjacent_scalar_indices_for_cross[k] = from_index[adjacent_vector_indices_h[6*j + k]];
            }
            if (j < NUMBER_OF_PENTAGONS)
                adjacent_scalar_indices_for_cross[5] = -1;
            else
            {
                if (to_index[adjacent_vector_indices_h[6*j + 5]] != j)
                    adjacent_scalar_indices_for_cross[5] = to_index[adjacent_vector_indices_h[6*j + 5]];
                else
                    adjacent_scalar_indices_for_cross[5] = from_index[adjacent_vector_indices_h[6*j + 5]];
            }
            retval = in_bool_calculator_long(adjacent_scalar_indices_for_cross, 6, from_index[i], &bool_0);
            retval = in_bool_calculator_long(adjacent_scalar_indices_for_cross, 6, to_index[i], &bool_1);
            if (bool_0 == 1 && bool_1 == 1)
            {
                if (first_found == 0)
                {
                    cell_0_for_cross = j;
                    first_found = 1;
                }
                else
                    cell_1_for_cross = j;
            }
        }
        counter = 0;
        for (int k = 0; k < 6; ++k)
        {
            if (to_index[adjacent_vector_indices_h[6*cell_0_for_cross + k]] == to_index[i] || from_index[adjacent_vector_indices_h[6*cell_0_for_cross + k]] == to_index[i] || to_index[adjacent_vector_indices_h[6*cell_0_for_cross + k]] == from_index[i] || from_index[adjacent_vector_indices_h[6*cell_0_for_cross + k]] == from_index[i])
            {
                face_of_cell_indices[counter] = k;
                ++counter;
            }
        }
        if (counter != 2)
            printf("Trouble detected, place 2.\n");
        recov_hor_par_pri_index[4*i] = adjacent_vector_indices_h[6*cell_0_for_cross + face_of_cell_indices[0]];
        recov_hor_par_pri_weight[4*i] = 2.0/3.0*0.5*(cos(direction[i] + M_PI/2)*cos(direction[recov_hor_par_pri_index[4*i]]) + sin(direction[i] + M_PI/2)*sin(direction[recov_hor_par_pri_index[4*i]]));
        recov_hor_par_pri_index[4*i + 1] = adjacent_vector_indices_h[6*cell_0_for_cross + face_of_cell_indices[1]];
        recov_hor_par_pri_weight[4*i + 1] = 2.0/3.0*0.5*(cos(direction[i] + M_PI/2)*cos(direction[recov_hor_par_pri_index[4*i + 1]]) + sin(direction[i] + M_PI/2)*sin(direction[recov_hor_par_pri_index[4*i + 1]]));
        counter = 0;
        for (int k = 0; k < 6; ++k)
        {
            if (to_index[adjacent_vector_indices_h[6*cell_1_for_cross + k]] == to_index[i] || from_index[adjacent_vector_indices_h[6*cell_1_for_cross + k]] == to_index[i] || to_index[adjacent_vector_indices_h[6*cell_1_for_cross + k]] == from_index[i] || from_index[adjacent_vector_indices_h[6*cell_1_for_cross + k]] == from_index[i])
            {
                face_of_cell_indices[counter] = k;
                ++counter;
            }
        }
        if (counter != 2)
            printf("Trouble detected, place 3.\n");
        recov_hor_par_pri_index[4*i + 2] = adjacent_vector_indices_h[6*cell_1_for_cross + face_of_cell_indices[0]];
        recov_hor_par_pri_weight[4*i + 2] = 2.0/3.0*0.5*(cos(direction[i] + M_PI/2)*cos(direction[recov_hor_par_pri_index[4*i + 2]]) + sin(direction[i] + M_PI/2)*sin(direction[recov_hor_par_pri_index[4*i + 2]]));
        recov_hor_par_pri_index[4*i + 3] = adjacent_vector_indices_h[6*cell_1_for_cross + face_of_cell_indices[1]];
        recov_hor_par_pri_weight[4*i + 3] = 2.0/3.0*0.5*(cos(direction[i] + M_PI/2)*cos(direction[recov_hor_par_pri_index[4*i + 3]]) + sin(direction[i] + M_PI/2)*sin(direction[recov_hor_par_pri_index[4*i + 3]]));
        sign = 1;
        find_angle_change(direction[i], direction_dual[i], &direction_change);
        if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
            sign = -1;
        for (int j = 0; j < 4; ++j)
        {
            if (j == 0)
            {
                h_curl_indices_dual[4*i + j] = i + NUMBER_OF_DUAL_VECTORS_PER_LAYER;
                h_curl_signs_dual[4*i + j] = sign;
            }
            if (j == 1)
            {
                if (sign == 1)
                    h_curl_indices_dual[4*i + j] = to_index_dual[i];
                else
                    h_curl_indices_dual[4*i + j] = from_index_dual[i];
                h_curl_signs_dual[4*i + j] = 1;
            }
            if (j == 2)
            {
                h_curl_indices_dual[4*i + j] = i;
                h_curl_signs_dual[4*i + j] = -sign;
            }
            if (j == 3)
            {
                if (sign == 1)
                    h_curl_indices_dual[4*i + j] = from_index_dual[i];
                else
                    h_curl_indices_dual[4*i + j] = to_index_dual[i];
                h_curl_signs_dual[4*i + j] = -1;
            }
            if (j == 0)
                recov_hor_ver_pri_index[4*i + j] = to_index[i];
            if (j == 1)
                recov_hor_ver_pri_index[4*i + j] = from_index[i];
            if (j == 2)
                recov_hor_ver_pri_index[4*i + j] = to_index[i] + NUMBER_OF_VECTORS_PER_LAYER;
            if (j == 3)
                recov_hor_ver_pri_index[4*i + j] = from_index[i] + NUMBER_OF_VECTORS_PER_LAYER;
            recov_hor_ver_pri_weight[4*i + j] = 0.25;
        }
    }
    free(face_of_cell_indices);
    free(adjacent_scalar_indices_for_cross);
    double weight_prefactor;
    for (int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
    {
        weight_prefactor = 2.0/6.0;
        if (i < NUMBER_OF_PENTAGONS)
            weight_prefactor = 2.0/5.0;
        for (int j = 0; j < 6; ++j)
        {
            recov_ver_0_pri_index[12*i + j] = adjacent_vector_indices_h[6*i + j];
            recov_ver_0_pri_weight[12*i + j] = 0.5*weight_prefactor*cos(direction[recov_ver_0_pri_index[12*i + j]]);
            recov_ver_0_dual_index[6*i + j] = adjacent_vector_indices_h[6*i + j];
            recov_ver_0_dual_weight[6*i + j] = weight_prefactor*cos(direction_dual[recov_ver_0_dual_index[6*i + j]]);
            recov_ver_1_pri_index[12*i + j] = adjacent_vector_indices_h[6*i + j];
            recov_ver_1_pri_weight[12*i + j] = 0.5*weight_prefactor*sin(direction[recov_ver_1_pri_index[12*i + j]]);
            recov_ver_1_dual_index[6*i + j] = adjacent_vector_indices_h[6*i + j];
            recov_ver_1_dual_weight[6*i + j] = weight_prefactor*sin(direction_dual[recov_ver_1_dual_index[6*i + j]]);
        }
        for (int j = 6; j < 12; ++j)
        {
            recov_ver_0_pri_index[12*i + j] = adjacent_vector_indices_h[6*i + j - 6] + NUMBER_OF_VECTORS_PER_LAYER;
            recov_ver_0_pri_weight[12*i + j] = 0.5*weight_prefactor*cos(direction[recov_ver_0_pri_index[12*i + j]]);
            recov_ver_1_pri_index[12*i + j] = adjacent_vector_indices_h[6*i + j - 6] + NUMBER_OF_VECTORS_PER_LAYER;
            recov_ver_1_pri_weight[12*i + j] = 0.5*weight_prefactor*sin(direction[recov_ver_1_pri_index[12*i + j]]);
        }
        if (i < NUMBER_OF_PENTAGONS)
        {
            recov_ver_0_pri_index[12*i + 5] = 0;
            recov_ver_0_pri_index[12*i + 11] = 0;
            recov_ver_1_pri_index[12*i + 5] = 0;
            recov_ver_1_pri_index[12*i + 11] = 0;
            recov_ver_0_dual_index[6*i + 5] = 0;
            recov_ver_1_dual_index[6*i + 5] = 0;
            recov_ver_0_pri_weight[12*i + 5] = 0;
            recov_ver_0_pri_weight[12*i + 11] = 0;
            recov_ver_1_pri_weight[12*i + 5] = 0;
            recov_ver_1_pri_weight[12*i + 11] = 0;
            recov_ver_0_dual_weight[6*i + 5] = 0;
            recov_ver_1_dual_weight[6*i + 5] = 0;
        }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_H; ++i)
    {
        sign = 1;
        find_angle_change(direction_dual[i], direction[i], &direction_change);
        if (rad2deg(direction_change) < -ORTH_CRITERION_DEG)
            sign = -1;
        for (int j = 0; j < 4; ++j)
        {
            if (j == 0)
            {
                h_curl_indices[4*i + j] = i + NUMBER_OF_VECTORS_PER_LAYER;
                h_curl_signs[4*i + j] = sign;
            }
            if (j == 1)
            {
                if (sign == 1)
                    h_curl_indices[4*i + j] = to_index[i];
                else
                    h_curl_indices[4*i + j] = from_index[i];
                h_curl_signs[4*i + j] = 1;
            }
            if (j == 2)
            {
                h_curl_indices[4*i + j] = i;
                h_curl_signs[4*i + j] = -sign;
            }
            if (j == 3)
            {
                if (sign == 1)
                    h_curl_indices[4*i + j] = from_index[i];
                else
                    h_curl_indices[4*i + j] = to_index[i];
                h_curl_signs[4*i + j] = -1;
            }
        }
    }
    free(direction_dual);
    int ncid_g_prop;
    if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid_g_prop)))
        ERR(retval);
    free(OUTPUT_FILE);
    int latitude_scalar_id, longitude_scalar_id, direction_id, latitude_vector_id, longitude_vector_id, latitude_scalar_dual_id, longitude_scalar_dual_id, z_scalar_id, z_vector_id, normal_distance_id, gravity_id, volume_id, area_id, recov_hor_par_dual_weight_id, recov_hor_ver_dual_weight_id, recov_hor_par_pri_weight_id, recov_hor_ver_pri_weight_id, recov_ver_0_pri_weight_id, recov_ver_0_dual_weight_id, recov_ver_1_pri_weight_id, recov_ver_1_dual_weight_id, z_vector_dual_id, normal_distance_dual_id, area_dual_id, f_vec_id, to_index_id, from_index_id, adjacent_vector_indices_h_id, vorticity_indices_id, h_curl_indices_id, recov_hor_par_dual_index_id, recov_hor_ver_dual_index_id, recov_hor_par_pri_index_id, recov_hor_ver_pri_index_id, recov_ver_0_pri_index_id, recov_ver_0_dual_index_id, recov_ver_1_pri_index_id, recov_ver_1_dual_index_id, to_index_dual_id, from_index_dual_id, vorticity_indices_dual_id, h_curl_indices_dual_id, adjacent_signs_h_id, vorticity_signs_id, h_curl_signs_id, vorticity_signs_dual_id, h_curl_signs_dual_id, vector_dual_one_layer_dimid;
    int scalar_dimid, scalar_h_dimid, scalar_dual_h_dimid, vector_dimid, scalar_h_dimid_6, vector_h_dimid, vector_h_dimid_11, vector_h_dimid_2, vector_h_dimid_4, vector_v_dimid_6, vector_dual_dimid, vector_dual_h_dimid, vector_dual_v_dimid_3, vector_v_dimid_12, vector_dual_h_dimid_4;
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_index", NUMBER_OF_SCALARS, &scalar_dimid)))
        ERR(retval);  
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_h_index", NUMBER_OF_SCALARS_H, &scalar_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_dual_h_index", NUMBER_OF_DUAL_SCALARS_H, &scalar_dual_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index", NUMBER_OF_VECTORS, &vector_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_index", NUMBER_OF_VECTORS_H, &vector_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "scalar_h_6_index", 6*NUMBER_OF_SCALARS_H, &scalar_h_dimid_6)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_11_index", 11*NUMBER_OF_VECTORS_H, &vector_h_dimid_11)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_2_index", 2*NUMBER_OF_VECTORS_H, &vector_h_dimid_2)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_h_4_index", 4*NUMBER_OF_VECTORS_H, &vector_h_dimid_4)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_v_6_index", 6*NUMBER_OF_VECTORS_V, &vector_v_dimid_6)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_v_12_index", 12*NUMBER_OF_VECTORS_V, &vector_v_dimid_12)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_index_h_dual", NUMBER_OF_DUAL_VECTORS_H, &vector_dual_h_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid_g_prop, "vector_dual_one_layer_dimid", NUMBER_OF_DUAL_VECTORS_PER_LAYER, &vector_dual_one_layer_dimid)))
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
    if ((retval = nc_def_var(ncid_g_prop, "z_scalar", NC_DOUBLE, 1, &scalar_dimid, &z_scalar_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_par_dual_weight", NC_DOUBLE, 1, &vector_h_dimid_2, &recov_hor_par_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_ver_dual_weight", NC_DOUBLE, 1, &vector_h_dimid_2, &recov_hor_ver_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_par_pri_weight", NC_DOUBLE, 1, &vector_h_dimid_4, &recov_hor_par_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_ver_pri_weight", NC_DOUBLE, 1, &vector_h_dimid_4, &recov_hor_ver_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_0_dual_weight", NC_DOUBLE, 1, &vector_v_dimid_6, &recov_ver_0_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_0_pri_weight", NC_DOUBLE, 1, &vector_v_dimid_12, &recov_ver_0_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_1_pri_weight", NC_DOUBLE, 1, &vector_v_dimid_12, &recov_ver_1_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_1_dual_weight", NC_DOUBLE, 1, &vector_v_dimid_6, &recov_ver_1_dual_weight_id)))
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
    if ((retval = nc_def_var(ncid_g_prop, "f_vec", NC_DOUBLE, 1, &vector_dual_one_layer_dimid, &f_vec_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid_g_prop, f_vec_id, "units", strlen("1/s"), "1/s")))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "direction", NC_DOUBLE, 1, &vector_h_dimid, &direction_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "latitude_vector", NC_DOUBLE, 1, &vector_h_dimid, &latitude_vector_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "longitude_vector", NC_DOUBLE, 1, &vector_h_dimid, &longitude_vector_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "to_index", NC_LONG, 1, &vector_h_dimid, &to_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "from_index", NC_LONG, 1, &vector_h_dimid, &from_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_vector_indices_h", NC_LONG, 1, &scalar_h_dimid_6, &adjacent_vector_indices_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_indices", NC_LONG, 1, &vector_dual_v_dimid_3, &vorticity_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "h_curl_indices", NC_LONG, 1, &vector_dual_h_dimid_4, &h_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_par_dual_index", NC_LONG, 1, &vector_h_dimid_2, &recov_hor_par_dual_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_ver_dual_index", NC_LONG, 1, &vector_h_dimid_2, &recov_hor_ver_dual_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_par_pri_index", NC_LONG, 1, &vector_h_dimid_4, &recov_hor_par_pri_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_hor_ver_pri_index", NC_LONG, 1, &vector_h_dimid_4, &recov_hor_ver_pri_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_0_pri_index", NC_LONG, 1, &vector_v_dimid_12, &recov_ver_0_pri_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_0_dual_index", NC_LONG, 1, &vector_v_dimid_6, &recov_ver_0_dual_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_1_pri_index", NC_LONG, 1, &vector_v_dimid_12, &recov_ver_1_pri_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "recov_ver_1_dual_index", NC_LONG, 1, &vector_v_dimid_6, &recov_ver_1_dual_index_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "to_index_dual", NC_LONG, 1, &vector_dual_h_dimid, &to_index_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "from_index_dual", NC_LONG, 1, &vector_dual_h_dimid, &from_index_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_indices_dual", NC_LONG, 1, &vector_v_dimid_6, &vorticity_indices_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "h_curl_indices_dual", NC_LONG, 1, &vector_h_dimid_4, &h_curl_indices_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "adjacent_signs_h", NC_SHORT, 1, &scalar_h_dimid_6, &adjacent_signs_h_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_signs", NC_SHORT, 1, &vector_dual_v_dimid_3, &vorticity_signs_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "h_curl_signs", NC_SHORT, 1, &vector_dual_h_dimid_4, &h_curl_signs_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "vorticity_signs_dual", NC_SHORT, 1, &vector_v_dimid_6, &vorticity_signs_dual_id)))
        ERR(retval);
    if ((retval = nc_def_var(ncid_g_prop, "h_curl_signs_dual", NC_SHORT, 1, &vector_h_dimid_4, &h_curl_signs_dual_id)))
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
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_0_pri_weight_id, &recov_ver_0_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_0_dual_weight_id, &recov_ver_0_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_1_pri_weight_id, &recov_ver_1_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid_g_prop, recov_ver_1_dual_weight_id, &recov_ver_1_dual_weight[0])))
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
    if ((retval = nc_put_var_long(ncid_g_prop, recov_ver_0_pri_index_id, &recov_ver_0_pri_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, recov_ver_0_dual_index_id, &recov_ver_0_dual_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, recov_ver_1_pri_index_id, &recov_ver_1_pri_index[0])))
        ERR(retval);
    if ((retval = nc_put_var_long(ncid_g_prop, recov_ver_1_dual_index_id, &recov_ver_1_dual_index[0])))
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
    if ((retval = nc_put_var_short(ncid_g_prop, vorticity_signs_dual_id, &vorticity_signs_dual[0])))
        ERR(retval);
    if ((retval = nc_put_var_short(ncid_g_prop, h_curl_signs_dual_id, &h_curl_signs_dual[0])))
        ERR(retval);
    if ((retval = nc_close(ncid_g_prop)))
        ERR(retval);
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
    free(recov_hor_par_dual_weight);
    free(recov_hor_ver_dual_weight);
    free(recov_hor_par_pri_weight);
    free(recov_hor_ver_pri_weight);
    free(recov_ver_0_pri_weight);
    free(recov_ver_0_dual_weight);
    free(recov_ver_1_pri_weight);
    free(recov_ver_1_dual_weight);
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
    free(recov_ver_0_pri_index);
    free(recov_ver_0_dual_index);
    free(recov_ver_1_pri_index);
    free(recov_ver_1_dual_index);
    free(to_index_dual);
    free(from_index_dual);
    free(vorticity_indices_dual);
    free(h_curl_indices_dual);
    free(adjacent_signs_h);
    free(vorticity_signs);
    free(h_curl_signs);
    free(vorticity_signs_dual);
    free(h_curl_signs_dual);
    free(area_dual);
    return 0;
}

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

int find_coords_from_triangle_on_face_index(long triangle_on_face_index, long res_id, long *coord_0, long *coord_1, long *coord_0_points_amount)
{
    short check = 1;
    long coord_1_pre = -1;
    long min_index, max_index, points_per_edge;
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

int find_triangle_on_face_index_from_coords(long coord_0, long coord_1, long res_id, long *triangle_on_face_index)
{
    int i = 0;
    *triangle_on_face_index = 0;
    long points_per_edge;
    int retval = find_points_per_edge(res_id, &points_per_edge);
    long coord_0_points_amount = points_per_edge;
    while (i < coord_1)
    {
        *triangle_on_face_index += coord_0_points_amount;
        coord_0_points_amount -= 1;
        ++i;
    }
    *triangle_on_face_index += coord_0;
    return 0;
}

int find_sigma_from_level(int level_index, double *result)
{
    *result = 1 - (level_index + 0.0)/NUMBER_OF_LAYERS;
    return 0;
}

int find_sigma_from_layer(int layer_index, double *result)
{
    *result = 1 - (layer_index + 0.5)/NUMBER_OF_LAYERS;
    return 0;
}


int find_triangle_indices_from_h_vector_index(int res_id, long i, long *point_0, long *point_1, long *point_2, long *point_3, long *point_4, long *point_5, long *dual_scalar_on_face_index, short *small_triangle_edge_index, short face_edges[][3], short face_vertices[][3], short edge_vertices [][2], short face_edges_reverse[][3])
{
    short face_index = (i - NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
    long on_face_index = i - (NUMBER_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
    int triangle_on_face_index = on_face_index/3;
    *small_triangle_edge_index = on_face_index - 3*triangle_on_face_index;
    int retval = find_triangle_edge_points(triangle_on_face_index, face_index, res_id, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, face_vertices, face_edges, face_edges_reverse);
    return retval;
}

int find_triangle_edge_points(long triangle_on_face_index, short face_index, long res_id, long *point_0, long *point_1, long *point_2, long *point_3, long *point_4, long *point_5, long *dual_scalar_on_face_index, short face_vertices[][3], short face_edges[][3], short face_edges_reverse[][3])
{
    long coord_0, coord_1, coord_0_points_amount;
    int retval = find_coords_from_triangle_on_face_index(triangle_on_face_index, res_id, &coord_0, &coord_1, &coord_0_points_amount);
    *dual_scalar_on_face_index = 1 + 2*triangle_on_face_index + coord_1;
    long points_per_edge, scalar_points_per_inner_face;
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

int find_triangle_on_face_index_from_dual_scalar_on_face_index(long dual_scalar_on_face_index, long res_id, long *triangle_on_face_index, short *points_downwards, short *special_case_bool, short *last_triangle_bool)
{
    short value_found = 0;
    int retval = 0;
    long triangle_on_face_index_pre, coord_0_pre, coord_1_pre, coord_0_points_amount_pre, dual_scalar_on_face_index_0, dual_scalar_on_face_index_1, dual_scalar_on_face_index_2, dual_scalar_on_face_index_3, points_per_edge;
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

int find_triangle_edge_points_from_dual_scalar_on_face_index(long dual_scalar_on_face_index, short face_index, long res_id, long *point_0, long *point_1, long *point_2, short face_vertices[][3], short face_edges[][3], short face_edges_reverse[][3])
{
    short points_downwards, special_case_bool, last_triangle_bool;
    long triangle_on_face_index, rhombuspoint_0, rhombuspoint_1, rhombuspoint_2, rhombuspoint_3, coord_0, coord_1, coord_0_points_amount, points_per_edge, dump, addpoint_0, addpoint_1;
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

int find_points_per_edge(long res_id, long *points_per_edge)
{
    *points_per_edge = (long) (pow(2, res_id) - 1);
    return 0;
}

int find_scalar_points_per_inner_face(long res_id, long *scalar_points_per_inner_face)
{
    *scalar_points_per_inner_face = (long) (0.5*(pow(2, res_id) - 2)*(pow(2, res_id) - 1));
    return 0;
}

int upscale_scalar_point(long res_id, long old_index, long *new_index)
{
    short edge_index, face_index;
    long points_per_edge, on_edge_index, scalar_points_per_inner_face, on_face_index, coord_0, coord_1, coord_0_points_amount;
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

int write_scalar_coordinates(long edgepoint_0, long edgepoint_1, long edgepoint_2, long point_0, long point_1, long point_2, short points_upwards, double x_unity[], double y_unity[], double z_unity[], double latitude_scalar[], double longitude_scalar[])
{
    double x_res, y_res, z_res, lat_res, lon_res;
    int retval = find_between_point(x_unity[edgepoint_0], y_unity[edgepoint_0], z_unity[edgepoint_0], x_unity[edgepoint_1], y_unity[edgepoint_1], z_unity[edgepoint_1], 0.5, &x_res, &y_res, &z_res);
    double test_0, test_1, test_2;
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

int find_triangles_per_face(long res_id, long *number_of_triangles_per_face)
{
    *number_of_triangles_per_face = (long) (pow(4, res_id));
    return 0;
}





