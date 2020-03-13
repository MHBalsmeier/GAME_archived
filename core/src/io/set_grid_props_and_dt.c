#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "../enum_and_typedefs.h"
#include "io.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int set_grid_properties(Grid *grid, Dualgrid *dualgrid, char GEO_PROP_FILE[])
{
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
    double *normal_distance_dual = malloc(NUMBER_OF_DUAL_VECTORS*sizeof(double));
    double *area_dual = malloc(NUMBER_OF_DUAL_VECTORS*sizeof(double));
    double *f_vec = malloc(NUMBER_OF_DUAL_VECTORS_PER_LAYER*sizeof(double));
    long *recov_ver_2_dual_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(long));
    long *recov_ver_1_dual_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(long));
    long *recov_ver_2_pri_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(long));
    long *recov_hor_ver_pri_index = malloc(4*NUMBER_OF_VECTORS_H*sizeof(long));
    long *recov_ver_1_pri_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(long));
    long *recov_hor_ver_dual_index = malloc(2*NUMBER_OF_VECTORS_H*sizeof(long));
    long *recov_hor_par_pri_index = malloc(2*NUMBER_OF_VECTORS_H*sizeof(long));
    long *to_index = malloc(NUMBER_OF_VECTORS_H*sizeof(long));
    long *from_index = malloc(NUMBER_OF_VECTORS_H*sizeof(long));
    long *adjacent_vector_indices_h = malloc(6*NUMBER_OF_SCALARS_H*sizeof(long));
    long *vorticity_indices = malloc(6*NUMBER_OF_SCALARS_H*sizeof(long));
    long *h_curl_indices = malloc(4*NUMBER_OF_VECTORS_H*sizeof(long));
    long *recov_hor_par_dual_index = malloc(11*NUMBER_OF_VECTORS_H*sizeof(long));
    long *to_index_dual = malloc(NUMBER_OF_DUAL_VECTORS_H*sizeof(long));
    long *from_index_dual = malloc(NUMBER_OF_DUAL_VECTORS_H*sizeof(long));
    long *vorticity_indices_dual = malloc(3*NUMBER_OF_DUAL_VECTORS_V*sizeof(long));
    long *h_curl_indices_dual = malloc(4*NUMBER_OF_DUAL_VECTORS_H*sizeof(long));
    short *adjacent_signs_h = (short *) malloc(6*NUMBER_OF_SCALARS_H*sizeof(short));
    short *vorticity_signs = (short *) malloc(6*NUMBER_OF_SCALARS_H*sizeof(short));
    short *h_curl_signs = (short *) malloc(4*NUMBER_OF_VECTORS_H*sizeof(short));
    short *vector_product_sign = (short *) malloc(NUMBER_OF_VECTORS_H*sizeof(short));
    short *vorticity_signs_dual = (short *) malloc(3*NUMBER_OF_DUAL_VECTORS_V*sizeof(short));
    short *h_curl_signs_dual = (short *) malloc(4*NUMBER_OF_DUAL_VECTORS_H*sizeof(short));
    int ncid;
    int retval;
    int normal_distance_id, gravity_id, volume_id, area_id, recov_hor_par_dual_weight_id, recov_hor_ver_dual_weight_id, recov_hor_par_pri_weight_id, recov_hor_ver_pri_weight_id, recov_ver_1_pri_weight_id, recov_ver_1_dual_weight_id, recov_ver_2_pri_weight_id, recov_ver_2_dual_weight_id, normal_distance_dual_id, area_dual_id, f_vec_id, to_index_id, from_index_id, adjacent_vector_indices_h_id, vorticity_indices_id, h_curl_indices_id, recov_hor_par_dual_index_id, recov_hor_ver_dual_index_id, recov_hor_par_pri_index_id, recov_hor_ver_pri_index_id, recov_ver_1_pri_index_id, recov_ver_1_dual_index_id, recov_ver_2_pri_index_id, recov_ver_2_dual_index_id, to_index_dual_id, from_index_dual_id, vorticity_indices_dual_id, h_curl_indices_dual_id, adjacent_signs_h_id, vorticity_signs_id, h_curl_signs_id, vector_product_sign_id, vorticity_signs_dual_id, h_curl_signs_dual_id;
    long vert_index, floor_index, h_index, layer_index;
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "normal_distance", &normal_distance_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "gravity", &gravity_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "volume", &volume_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "area", &area_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_hor_par_dual_weight", &recov_hor_par_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_hor_ver_dual_weight", &recov_hor_ver_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_hor_par_pri_weight", &recov_hor_par_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_hor_ver_pri_weight", &recov_hor_ver_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_1_dual_weight", &recov_ver_1_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_1_pri_weight", &recov_ver_1_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_2_pri_weight", &recov_ver_2_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_2_dual_weight", &recov_ver_2_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "normal_distance_dual", &normal_distance_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "area_dual", &area_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "f_vec", &f_vec_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "to_index", &to_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "from_index", &from_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_vector_indices_h", &adjacent_vector_indices_h_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vorticity_indices", &vorticity_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "h_curl_indices", &h_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_hor_par_dual_index", &recov_hor_par_dual_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_hor_ver_dual_index", &recov_hor_ver_dual_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_hor_par_pri_index", &recov_hor_par_pri_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_hor_ver_pri_index", &recov_hor_ver_pri_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_1_pri_index", &recov_ver_1_pri_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_1_dual_index", &recov_ver_1_dual_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_2_pri_index", &recov_ver_2_pri_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_2_dual_index", &recov_ver_2_dual_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "to_index_dual", &to_index_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "from_index_dual", &from_index_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vorticity_indices_dual", &vorticity_indices_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "h_curl_indices_dual", &h_curl_indices_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_signs_h", &adjacent_signs_h_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vorticity_signs", &vorticity_signs_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "h_curl_signs", &h_curl_signs_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vector_product_sign", &vector_product_sign_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vorticity_signs_dual", &vorticity_signs_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "h_curl_signs_dual", &h_curl_signs_dual_id)))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, normal_distance_id, &normal_distance[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, gravity_id, &gravity[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, volume_id, &volume[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, area_id, &area[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_hor_par_dual_weight_id, &recov_hor_par_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_hor_ver_dual_weight_id, &recov_hor_ver_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_hor_par_pri_weight_id, &recov_hor_par_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_hor_ver_pri_weight_id, &recov_hor_ver_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_ver_1_pri_weight_id, &recov_ver_1_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_ver_1_dual_weight_id, &recov_ver_1_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_ver_2_pri_weight_id, &recov_ver_2_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_ver_2_dual_weight_id, &recov_ver_2_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, normal_distance_dual_id, &normal_distance_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, area_dual_id, &area_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, f_vec_id, &f_vec[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, to_index_id, &to_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, from_index_id, &from_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_vector_indices_h_id, &adjacent_vector_indices_h[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, vorticity_indices_id, &vorticity_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, h_curl_indices_id, &h_curl_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, recov_hor_par_dual_index_id, &recov_hor_par_dual_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, recov_hor_ver_dual_index_id, &recov_hor_ver_dual_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, recov_hor_par_pri_index_id, &recov_hor_par_pri_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, recov_hor_ver_pri_index_id, &recov_hor_ver_pri_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, recov_ver_1_pri_index_id, &recov_ver_1_pri_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, recov_ver_1_dual_index_id, &recov_ver_1_dual_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, recov_ver_2_pri_index_id, &recov_ver_2_pri_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, recov_ver_2_dual_index_id, &recov_ver_2_dual_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, to_index_dual_id, &to_index_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, from_index_dual_id, &from_index_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, vorticity_indices_dual_id, &vorticity_indices_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, h_curl_indices_dual_id, &h_curl_indices_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_short(ncid, adjacent_signs_h_id, &adjacent_signs_h[0])))
        ERR(retval);
    if ((retval = nc_get_var_short(ncid, vorticity_signs_id, &vorticity_signs[0])))
        ERR(retval);
    if ((retval = nc_get_var_short(ncid, h_curl_signs_id, &h_curl_signs[0])))
        ERR(retval);
    if ((retval = nc_get_var_short(ncid, vector_product_sign_id, &vector_product_sign[0])))
        ERR(retval);
    if ((retval = nc_get_var_short(ncid, vorticity_signs_dual_id, &vorticity_signs_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_short(ncid, h_curl_signs_dual_id, &h_curl_signs_dual[0])))
        ERR(retval);
    if ((retval = nc_close(ncid)))
        ERR(retval);
    int i, j, k;
    for (i = 0; i < NUMBER_OF_SCALARS_H; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            grid -> adjacent_vector_indices_h[6*i + j] = adjacent_vector_indices_h[6*i + j];
            grid -> adjacent_signs_h[6*i + j] = adjacent_signs_h[6*i + j];
            grid -> vorticity_indices[6*i + j] = vorticity_indices[6*i +j ];
            grid -> vorticity_signs[6*i + j] = vorticity_signs[6*i + j];
        }
    }
    for (i = 0; i < NUMBER_OF_SCALARS; ++i)
        grid -> volume[i] = volume[i];
    for (i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        grid -> to_index[i] = to_index[i];
        grid -> from_index[i] = from_index[i];
        for (int j = 0; j < 4; ++j)
        {
            grid -> h_curl_indices[4*i + j] = h_curl_indices[4*i + j];
            grid -> h_curl_signs[4*i + j] = h_curl_signs[4*i + j];
            grid -> adjacent_signs_h[4*i + j] = adjacent_signs_h[4*i + j];
            grid -> vorticity_indices[4*i + j] = vorticity_indices[4*i + j];
        }
    }
    for (i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        grid -> gravity[i] = gravity[i];
        grid -> normal_distance[i] = normal_distance[i];
        grid -> area[i] = area[i];
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        grid -> vector_product_sign[i] = vector_product_sign[i];
        for (int j = 0; j < 11; ++j)
        {
            grid -> recov_hor_par_dual_index[11*i + j] = recov_hor_par_dual_index[11*i + j];
            grid -> recov_hor_par_dual_weight[11*i + j] = recov_hor_par_dual_weight[11*i + j];
        }
        for (int j = 0; j < 2; ++j)
        {
            grid -> recov_hor_ver_dual_index[2*i + j] = recov_hor_ver_dual_index[2*i + j];
            grid -> recov_hor_ver_dual_weight[2*i + j] = recov_hor_ver_dual_weight[2*i + j];
            grid -> recov_hor_par_pri_index[2*i + j] = recov_hor_par_pri_index[2*i + j];
            grid -> recov_hor_par_pri_weight[2*i + j] = recov_hor_par_pri_weight[2*i + j];
        }
        for (int j = 0; j < 4; ++j)
        {
            grid -> recov_hor_ver_pri_index[4*i + j] = recov_hor_ver_pri_index[4*i + j];
            grid -> recov_hor_ver_pri_weight[4*i + j] = recov_hor_ver_pri_weight[4*i + j];
        }
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            grid -> recov_ver_1_pri_index[6*i + j] = recov_ver_1_pri_index[6*i + j];
            grid -> recov_ver_1_pri_weight[6*i + j] = recov_ver_1_pri_weight[6*i + j];
            grid -> recov_ver_1_dual_index[6*i + j] = recov_ver_1_dual_index[6*i + j];
            grid -> recov_ver_1_dual_weight[6*i + j] = recov_ver_1_dual_weight[6*i + j];
            grid -> recov_ver_2_pri_index[6*i + j] = recov_ver_2_pri_index[6*i + j];
            grid -> recov_ver_2_pri_weight[6*i + j] = recov_ver_2_pri_weight[6*i + j];
            grid -> recov_ver_2_dual_index[6*i + j] = recov_ver_2_dual_index[6*i + j];
            grid -> recov_ver_2_dual_weight[6*i + j] = recov_ver_2_dual_weight[6*i + j];
        }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_PER_LAYER; ++i)
        dualgrid -> f_vec[i] = f_vec[i];
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; ++i)
    {
        dualgrid -> area[i] = area_dual[i];
        dualgrid -> normal_distance[i] = normal_distance_dual[i];
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_H; ++i)
    {
        dualgrid -> to_index[i] = to_index_dual[i];
        dualgrid -> from_index[i] = from_index_dual[i];
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_H; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            dualgrid -> h_curl_indices[4*i + j] = h_curl_indices_dual[4*i + j];
            dualgrid -> h_curl_signs[4*i + j] = h_curl_signs_dual[4*i + j];
        }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_V; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            dualgrid -> vorticity_indices[3*i + j] = vorticity_indices_dual[3*i + j];
            dualgrid -> vorticity_signs[3*i + j] = vorticity_signs_dual[3*i + j];
        }
    }
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
    free(normal_distance_dual);
    free(area_dual);
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
    return 0;
}

double calc_delta_t(int res_id)
{
    double delta_t;
    if (res_id == 2)
        delta_t = 1800;
    if (res_id == 3)
        delta_t = 900;
    if (res_id == 4)
        delta_t = 450;
    if (res_id == 5)
        delta_t = 225;
    return delta_t;
}
