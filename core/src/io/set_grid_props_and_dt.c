#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "../enum_and_typedefs.h"
#include "io.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

Grid grid;
Dualgrid dualgrid;

int set_grid_properties(char GEO_PROP_FILE[])
{
    double *test;
    test = (double *) calloc(1, sizeof(double));
    double *latitude_scalar, *longitude_scalar, *z_scalar, *latitude_vector, *longitude_vector, *z_vector, *parallel_distance, *normal_distance, *gravity, *volume, *area, *recov_hor_par_dual_weight, *recov_hor_ver_dual_weight, *recov_hor_par_pri_weight, *recov_hor_ver_pri_weight, *recov_ver_1_pri_weight, *recov_ver_1_dual_weight, *recov_ver_2_pri_weight, *recov_ver_2_dual_weight, *latitude_scalar_dual, *longitude_scalar_dual, *z_scalar_dual, *latitude_vector_dual, *longitude_vector_dual, *z_vector_dual, *normal_distance_dual, *area_dual, *f_vec, *parallel_distance_dual;
    long *to_indices, *from_indices, *adjacent_vector_indices_h, *adjacent_vector_index_lower, *adjacent_vector_index_upper, *adjacent_scalar_indices_h, *adjacent_scalar_index_lower, *adjacent_scalar_index_upper, *vorticity_indices, *h_curl_indices, *recov_hor_par_dual_index, *recov_hor_ver_dual_index, *recov_hor_par_pri_index, *recov_hor_ver_pri_index, *recov_ver_1_pri_index, *recov_ver_1_dual_index, *recov_ver_2_pri_index, *recov_ver_2_dual_index, *to_indices_dual, *from_indices_dual, *adjacent_vector_indices_h_dual, *adjacent_vector_index_upper_dual, *adjacent_vector_index_lower_dual, *adjacent_scalar_indices_h_dual, *adjacent_scalar_index_lower_dual, *adjacent_scalar_index_upper_dual, *vorticity_indices_dual, *h_curl_indices_dual;
    short  *adjacent_signs_h, *vorticity_signs, *h_curl_signs, *vector_product_sign, *adjacent_signs_h_dual, *vorticity_signs_dual, *h_curl_signs_dual;   latitude_scalar = (double *) calloc(NUMBER_OF_SCALARS, sizeof(double));
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
    int ncid;
    int retval;
    int latitude_scalar_id, longitude_scalar_id, z_scalar_id, latitude_vector_id, longitude_vector_id, z_vector_id, parallel_distance_id, normal_distance_id, gravity_id, volume_id, area_id, recov_hor_par_dual_weight_id, recov_hor_ver_dual_weight_id, recov_hor_par_pri_weight_id, recov_hor_ver_pri_weight_id, recov_ver_1_pri_weight_id, recov_ver_1_dual_weight_id, recov_ver_2_pri_weight_id, recov_ver_2_dual_weight_id, latitude_scalar_dual_id, longitude_scalar_dual_id, z_scalar_dual_id, latitude_vector_dual_id, longitude_vector_dual_id, z_vector_dual_id, normal_distance_dual_id, area_dual_id, f_vec_id, parallel_distance_dual_id, to_indices_id, from_indices_id, adjacent_vector_indices_h_id, adjacent_vector_index_lower_id, adjacent_vector_index_upper_id, adjacent_scalar_indices_h_id, adjacent_scalar_index_lower_id, adjacent_scalar_index_upper_id, vorticity_indices_id, h_curl_indices_id, recov_hor_par_dual_index_id, recov_hor_ver_dual_index_id, recov_hor_par_pri_index_id, recov_hor_ver_pri_index_id, recov_ver_1_pri_index_id, recov_ver_1_dual_index_id, recov_ver_2_pri_index_id, recov_ver_2_dual_index_id, to_indices_dual_id, from_indices_dual_id, adjacent_vector_indices_h_dual_id, adjacent_vector_index_upper_dual_id, adjacent_vector_index_lower_dual_id, adjacent_scalar_indices_h_dual_id, adjacent_scalar_index_lower_dual_id, adjacent_scalar_index_upper_dual_id, vorticity_indices_dual_id, h_curl_indices_dual_id, adjacent_signs_h_id, vorticity_signs_id, h_curl_signs_id, vector_product_sign_id, adjacent_signs_h_dual_id, vorticity_signs_dual_id, h_curl_signs_dual_id;
    long vert_index, floor_index, h_index, layer_index;
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "latitude_scalar", &latitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_scalar", &longitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_scalar", &z_scalar_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "latitude_vector", &latitude_vector_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_vector", &longitude_vector_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_vector", &z_vector_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "parallel_distance", &parallel_distance_id)))
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
    if ((retval = nc_inq_varid(ncid, "latitude_scalar_dual", &latitude_scalar_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_scalar_dual", &longitude_scalar_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_scalar_dual", &z_scalar_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "latitude_vector_dual", &latitude_vector_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_vector_dual", &longitude_vector_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_vector_dual", &z_vector_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "normal_distance_dual", &normal_distance_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "area_dual", &area_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "f_vec", &f_vec_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "parallel_distance_dual", &parallel_distance_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "to_indices", &to_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "from_indices", &from_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_vector_indices_h", &adjacent_vector_indices_h_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_vector_index_lower", &adjacent_vector_index_lower_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_vector_index_upper", &adjacent_vector_index_upper_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_scalar_indices_h", &adjacent_scalar_indices_h_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_scalar_index_lower", &adjacent_scalar_index_lower_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_scalar_index_upper", &adjacent_scalar_index_upper_id)))
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
    if ((retval = nc_inq_varid(ncid, "to_indices_dual", &to_indices_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "from_indices_dual", &from_indices_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_vector_indices_h_dual", &adjacent_vector_indices_h_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_vector_index_upper_dual", &adjacent_vector_index_upper_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_vector_index_lower_dual", &adjacent_vector_index_lower_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_scalar_indices_h_dual", &adjacent_scalar_indices_h_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_scalar_index_lower_dual", &adjacent_scalar_index_lower_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_scalar_index_upper_dual", &adjacent_scalar_index_upper_dual_id)))
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
    if ((retval = nc_inq_varid(ncid, "adjacent_signs_h_dual", &adjacent_signs_h_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vorticity_signs_dual", &vorticity_signs_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "h_curl_signs_dual", &h_curl_signs_dual_id)))
        ERR(retval);
     if ((retval = nc_get_var_double(ncid, latitude_scalar_id, &latitude_scalar[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_scalar_id, &longitude_scalar[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_scalar_id, &z_scalar[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_vector_id, &latitude_vector[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_vector_id, &longitude_vector[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_vector_id, &z_vector[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, parallel_distance_id, &parallel_distance[0])))
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
    if ((retval = nc_get_var_double(ncid, latitude_scalar_dual_id, &latitude_scalar_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_scalar_dual_id, &longitude_scalar_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_scalar_dual_id, &z_scalar_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_vector_dual_id, &latitude_vector_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_vector_dual_id, &longitude_vector_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_vector_dual_id, &z_vector_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, normal_distance_dual_id, &normal_distance_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, area_dual_id, &area_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, f_vec_id, &f_vec[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, parallel_distance_dual_id, &parallel_distance_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, to_indices_id, &to_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, from_indices_id, &from_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_vector_indices_h_id, &adjacent_vector_indices_h[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_vector_index_lower_id, &adjacent_vector_index_lower[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_vector_index_upper_id, &adjacent_vector_index_upper[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_scalar_indices_h_id, &adjacent_scalar_indices_h[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_scalar_index_lower_id, &adjacent_scalar_index_lower[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_scalar_index_upper_id, &adjacent_scalar_index_upper[0])))
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
    if ((retval = nc_get_var_long(ncid, to_indices_dual_id, &to_indices_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, from_indices_dual_id, &from_indices_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_vector_indices_h_dual_id, &adjacent_vector_indices_h_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_vector_index_upper_dual_id, &adjacent_vector_index_upper_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_vector_index_lower_dual_id, &adjacent_vector_index_lower_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_scalar_indices_h_dual_id, &adjacent_scalar_indices_h_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_scalar_index_lower_dual_id, &adjacent_scalar_index_lower_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(ncid, adjacent_scalar_index_upper_dual_id, &adjacent_scalar_index_upper_dual[0])))
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
    if ((retval = nc_get_var_short(ncid, adjacent_signs_h_dual_id, &adjacent_signs_h_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_short(ncid, vorticity_signs_dual_id, &vorticity_signs_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_short(ncid, h_curl_signs_dual_id, &h_curl_signs_dual[0])))
        ERR(retval);
    if ((retval = nc_close(ncid)))
        ERR(retval);
    int i, j, k;
    for (i = 0; i < NUMBER_OF_SCALARS; i++)
    {
        grid.latitude_scalar[i] = latitude_scalar[i];
        grid.longitude_scalar[i] = longitude_scalar[i];
        grid.z_scalar[i] = z_scalar[i];
        grid.volume[i] = volume[i];
        grid.adjacent_scalar_index_lower[i] = adjacent_scalar_index_lower[i];
        grid.adjacent_scalar_index_upper[i] = adjacent_scalar_index_upper[i];
        grid.adjacent_vector_index_lower[i] = adjacent_vector_index_lower[i];
        grid.adjacent_vector_index_upper[i] = adjacent_vector_index_upper[i];
        for (int j = 0; j < 6; j++)
        {
            grid.adjacent_scalar_indices_h[6*i + j] = adjacent_scalar_indices_h[6*i + j];
            grid.adjacent_vector_indices_h[6*i + j] = adjacent_vector_indices_h[6*i + j];
            grid.adjacent_signs_h[6*i + j] = adjacent_signs_h[6*i + j];
            grid.vorticity_indices[6*i + j] = vorticity_indices[6*i +j ];
            grid.vorticity_signs[6*i + j] = vorticity_signs[6*i + j];
        }
    }
    for (i = 0; i < NUMBER_OF_VECTORS; i++)
    {		
        vert_index = i/(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H);
		floor_index = vert_index*(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H);
		h_index = i - floor_index;
        grid.latitude_vector[i] = latitude_vector[i];
        grid.longitude_vector[i] = longitude_vector[i];
        grid.z_vector[i] = z_vector[i];
        grid.normal_distance[i] = normal_distance[i];
        grid.area[i] = area[i];
        grid.to_indices[i] = to_indices[i];
        grid.from_indices[i] = from_indices[i];
        for (int j = 0; j < 2; j++)
           grid.parallel_distance[2*i + j] = parallel_distance[2*i + j];
        for (int j = 0; j < 4; j++)
        {
            grid.h_curl_indices[4*i + j] = h_curl_indices[4*i + j];
            grid.h_curl_signs[4*i + j] = h_curl_signs[4*i + j];
            grid.adjacent_signs_h[4*i + j] = adjacent_signs_h[4*i + j];
            grid.vorticity_indices[4*i + j] = vorticity_indices[4*i + j];
        }
    }
    for (int i = 0; i < NUMBER_OF_H_VECTORS; i++)
    {
        grid.vector_product_sign[i] = vector_product_sign[i];
        for (int j = 0; j < 11; j++)
        {
            grid.recov_hor_par_dual_index[11*i + j] = recov_hor_par_dual_index[11*i + j];
            grid.recov_hor_par_dual_weight[11*i + j] = recov_hor_par_dual_weight[11*i + j];
        }
        for (int j = 0; j < 2; j++)
        {
            grid.recov_hor_ver_dual_index[2*i + j] = recov_hor_ver_dual_index[2*i + j];
            grid.recov_hor_ver_dual_weight[2*i + j] = recov_hor_ver_dual_weight[2*i + j];
            grid.recov_hor_par_pri_index[2*i + j] = recov_hor_par_pri_index[2*i + j];
            grid.recov_hor_par_pri_weight[2*i + j] = recov_hor_par_pri_weight[2*i + j];
        }
        for (int j = 0; j < 4; j++)
        {
            grid.recov_hor_ver_pri_index[4*i + j] = recov_hor_ver_pri_index[4*i + j];
            grid.recov_hor_ver_pri_weight[4*i + j] = recov_hor_ver_pri_weight[4*i + j];
        }
    }
    for (int i = 0; i < NUMBER_OF_V_VECTORS; i++)
    {
        grid.gravity[i] = gravity[i];
        for (int j = 0; j < 6; j++)
        {
            grid.recov_ver_1_pri_index[6*i + j] = recov_ver_1_pri_index[6*i + j];
            grid.recov_ver_1_pri_weight[6*i + j] = recov_ver_1_pri_weight[6*i + j];
            grid.recov_ver_1_dual_index[6*i + j] = recov_ver_1_dual_index[6*i + j];
            grid.recov_ver_1_dual_weight[6*i + j] = recov_ver_1_dual_weight[6*i + j];
            grid.recov_ver_2_pri_index[6*i + j] = recov_ver_2_pri_index[6*i + j];
            grid.recov_ver_2_pri_weight[6*i + j] = recov_ver_2_pri_weight[6*i + j];
            grid.recov_ver_2_dual_index[6*i + j] = recov_ver_2_dual_index[6*i + j];
            grid.recov_ver_2_dual_weight[6*i + j] = recov_ver_2_dual_weight[6*i + j];
        }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_SCALARS; i++)
    {
        dualgrid.latitude_scalar[i] = latitude_scalar_dual[i];
        dualgrid.longitude_scalar[i] = longitude_scalar_dual[i];
        dualgrid.z_scalar[i] = z_scalar_dual[i];
        dualgrid.adjacent_scalar_index_lower[i] = adjacent_scalar_index_lower_dual[i];
        dualgrid.adjacent_scalar_index_upper[i] = adjacent_scalar_index_upper_dual[i];
        dualgrid.adjacent_vector_index_upper[i] = adjacent_vector_index_upper_dual[i];
        dualgrid.adjacent_vector_index_lower[i] = adjacent_vector_index_lower_dual[i];
            for (int k = 0; k < 3; k++)
            {
                dualgrid.adjacent_scalar_indices_h[3*i + k] = adjacent_scalar_indices_h_dual[3*i + k];
                dualgrid.adjacent_vector_indices_h[3*i + k] = adjacent_vector_indices_h_dual[3*i + k];
                dualgrid.adjacent_signs_h[3*i + k] = adjacent_signs_h_dual[3*i + k];
            }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; i++)
    {
        dualgrid.latitude_vector[i] = latitude_vector_dual[i];
        dualgrid.longitude_vector[i] = longitude_vector_dual[i];
        dualgrid.z_vector[i] = z_vector_dual[i];
        dualgrid.area[i] = area_dual[i];
        dualgrid.f_vec[i] = f_vec[i];
        dualgrid.normal_distance[i] = normal_distance_dual[i];
        dualgrid.to_indices[i] = to_indices_dual[i];
        dualgrid.from_indices[i] = from_indices_dual[i];
        for (int j = 0; j < 2; j++)
            dualgrid.parallel_distance[2*i + j] = parallel_distance_dual[2*i + j];
    }
    for (int i = 0; i < NUMBER_OF_DUAL_H_VECTORS; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            dualgrid.h_curl_indices[4*i + j] = h_curl_indices_dual[4*i + j];
            dualgrid.h_curl_signs[4*i + j] = h_curl_signs_dual[4*i + j];
        }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_V_VECTORS; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            dualgrid.vorticity_indices[3*i + j] = vorticity_indices_dual[3*i + j];
            dualgrid.vorticity_signs[3*i + j] = vorticity_signs_dual[3*i + j];
        }
    }
    // freeing the memory 
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

double calc_delta_t(int res_id)
{
    double delta_t;
    if (res_id == 2)
        delta_t = 900;
    return delta_t;
}
