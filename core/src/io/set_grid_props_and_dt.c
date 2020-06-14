#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "../enum_and_typedefs.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

void grid_check_failed();

int set_grid_properties(Grid *grid, Dualgrid *dualgrid, char GEO_PROP_FILE[])
{
    double *normal_distance = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *volume = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *area = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *z_scalar = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *z_vector = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *recov_hor_par_dual_weight = malloc(2*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_hor_ver_dual_weight = malloc(2*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_hor_par_pri_weight = malloc(10*NUMBER_OF_VECTORS_H*sizeof(double));
    double *recov_ver_0_pri_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *recov_ver_1_pri_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *recov_ver_0_dual_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *recov_ver_1_dual_weight = malloc(6*NUMBER_OF_VECTORS_V*sizeof(double));
    double *normal_distance_dual = malloc(NUMBER_OF_DUAL_VECTORS*sizeof(double));
    double *area_dual = malloc((NUMBER_OF_DUAL_H_VECTORS + NUMBER_OF_H_VECTORS)*sizeof(double));
    double *f_vec = malloc(NUMBER_OF_DUAL_VECTORS_PER_LAYER*sizeof(double));
    double *direction = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *gravity_potential = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *gravity = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *vertical_contravar_unit = malloc(3*NUMBER_OF_VECTORS_V*(NUMBER_OF_ORO_LAYERS + 1)*sizeof(double));
    int *to_index = malloc(NUMBER_OF_VECTORS_H*sizeof(int));
    int *from_index = malloc(NUMBER_OF_VECTORS_H*sizeof(int));
    int *recov_ver_0_pri_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(int));
    int *recov_ver_1_pri_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(int));
    int *recov_ver_0_dual_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(int));
    int *recov_ver_1_dual_index = malloc(6*NUMBER_OF_VECTORS_V*sizeof(int));
    int *recov_hor_ver_pri_index = malloc(4*NUMBER_OF_VECTORS_H*sizeof(int));
    int *recov_hor_ver_dual_index = malloc(2*NUMBER_OF_VECTORS_H*sizeof(int));
    int *recov_hor_par_pri_index = malloc(10*NUMBER_OF_VECTORS_H*sizeof(int));
    int *adjacent_vector_indices_h = malloc(6*NUMBER_OF_SCALARS_H*sizeof(int));
    int *vorticity_indices = malloc(4*NUMBER_OF_VECTORS_H*sizeof(int));
    int *h_curl_indices = malloc(4*NUMBER_OF_DUAL_VECTORS_H*sizeof(int));
    int *recov_hor_par_dual_index = malloc(2*NUMBER_OF_VECTORS_H*sizeof(int));
    int *adjacent_signs_h = malloc(6*NUMBER_OF_SCALARS_H*sizeof(int));
    int *vorticity_signs = malloc(4*NUMBER_OF_VECTORS_H*sizeof(int));
    int *h_curl_signs = malloc(4*NUMBER_OF_DUAL_VECTORS_H*sizeof(int));
    int *adjacent_vector_indices_dual_h = malloc(3*NUMBER_OF_DUAL_SCALARS_H*sizeof(int));
    int ncid, retval;
    int normal_distance_id, volume_id, area_id, z_scalar_id, z_vector_id, recov_hor_par_dual_weight_id, recov_hor_ver_dual_weight_id, recov_hor_par_pri_weight_id, recov_ver_0_pri_weight_id, recov_ver_0_dual_weight_id, recov_ver_1_pri_weight_id, recov_ver_1_dual_weight_id, normal_distance_dual_id, area_dual_id, f_vec_id, to_index_id, from_index_id, adjacent_vector_indices_h_id, vorticity_indices_id, h_curl_indices_id, recov_hor_par_dual_index_id, recov_hor_ver_dual_index_id, recov_hor_par_pri_index_id, recov_hor_ver_pri_index_id, recov_ver_0_pri_index_id, recov_ver_0_dual_index_id, recov_ver_1_pri_index_id, recov_ver_1_dual_index_id, adjacent_signs_h_id, vorticity_signs_id, h_curl_signs_id, direction_id, gravity_potential_id, gravity_id, vertical_contravar_unit_id, adjacent_vector_indices_dual_h_id;
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "normal_distance", &normal_distance_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "volume", &volume_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "area", &area_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_scalar", &z_scalar_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "gravity_potential", &gravity_potential_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "gravity", &gravity_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vertical_contravar_unit", &vertical_contravar_unit_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_vector", &z_vector_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_hor_par_dual_weight", &recov_hor_par_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_hor_ver_dual_weight", &recov_hor_ver_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_hor_par_pri_weight", &recov_hor_par_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_0_dual_weight", &recov_ver_0_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_0_pri_weight", &recov_ver_0_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_1_pri_weight", &recov_ver_1_pri_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_1_dual_weight", &recov_ver_1_dual_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "normal_distance_dual", &normal_distance_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "area_dual", &area_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "f_vec", &f_vec_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "to_index", &to_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "direction", &direction_id)))
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
    if ((retval = nc_inq_varid(ncid, "recov_ver_0_pri_index", &recov_ver_0_pri_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_0_dual_index", &recov_ver_0_dual_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_1_pri_index", &recov_ver_1_pri_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_1_dual_index", &recov_ver_1_dual_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_signs_h", &adjacent_signs_h_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vorticity_signs", &vorticity_signs_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "h_curl_signs", &h_curl_signs_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_vector_indices_dual_h", &adjacent_vector_indices_dual_h_id)))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, normal_distance_id, &normal_distance[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, volume_id, &volume[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, area_id, &area[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_scalar_id, &z_scalar[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, gravity_potential_id, &gravity_potential[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, gravity_id, &gravity[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, vertical_contravar_unit_id, &vertical_contravar_unit[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_vector_id, &z_vector[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_hor_par_dual_weight_id, &recov_hor_par_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_hor_ver_dual_weight_id, &recov_hor_ver_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_hor_par_pri_weight_id, &recov_hor_par_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_ver_0_pri_weight_id, &recov_ver_0_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_ver_0_dual_weight_id, &recov_ver_0_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_ver_1_pri_weight_id, &recov_ver_1_pri_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_ver_1_dual_weight_id, &recov_ver_1_dual_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, normal_distance_dual_id, &normal_distance_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, area_dual_id, &area_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, direction_id, &direction[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, f_vec_id, &f_vec[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, to_index_id, &to_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, from_index_id, &from_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, adjacent_vector_indices_h_id, &adjacent_vector_indices_h[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, adjacent_vector_indices_dual_h_id, &adjacent_vector_indices_dual_h[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, vorticity_indices_id, &vorticity_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, h_curl_indices_id, &h_curl_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, recov_hor_par_dual_index_id, &recov_hor_par_dual_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, recov_hor_ver_dual_index_id, &recov_hor_ver_dual_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, recov_hor_par_pri_index_id, &recov_hor_par_pri_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, recov_hor_ver_pri_index_id, &recov_hor_ver_pri_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, recov_ver_0_pri_index_id, &recov_ver_0_pri_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, recov_ver_0_dual_index_id, &recov_ver_0_dual_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, recov_ver_1_pri_index_id, &recov_ver_1_pri_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, recov_ver_1_dual_index_id, &recov_ver_1_dual_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, adjacent_signs_h_id, &adjacent_signs_h[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, vorticity_signs_id, &vorticity_signs[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, h_curl_signs_id, &h_curl_signs[0])))
        ERR(retval);
    if ((retval = nc_close(ncid)))
        ERR(retval);
    for (int i = 0; i < NUMBER_OF_SCALARS_H; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            grid -> adjacent_vector_indices_h[6*i + j] = adjacent_vector_indices_h[6*i + j];
            if (grid -> adjacent_vector_indices_h[6*i + j] >= NUMBER_OF_VECTORS_H || grid -> adjacent_vector_indices_h[6*i + j] < -1)
                grid_check_failed();
            grid -> adjacent_signs_h[6*i + j] = adjacent_signs_h[6*i + j];
            if (grid -> adjacent_signs_h[6*i + j] != -1 && grid -> adjacent_signs_h[6*i + j] != 1)
            {
            	if (i >= NUMBER_OF_PENTAGONS || j != 5 || grid -> adjacent_signs_h[6*i + j] != 0)
                	grid_check_failed();
        	}
        }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_SCALARS_H; ++i)
    {
    	for (int j = 0; j < 3; ++j)
    	{
    	dualgrid -> adjacent_vector_indices_h[3*i + j] = adjacent_vector_indices_dual_h[3*i + j];
    	if (dualgrid -> adjacent_vector_indices_h[3*i + j] < 0 || dualgrid -> adjacent_vector_indices_h[3*i + j] >= NUMBER_OF_VECTORS_H)
    		grid_check_failed();
    	}
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            dualgrid -> vorticity_indices[4*i + j] = vorticity_indices[4*i +j];
            if (dualgrid -> vorticity_indices[4*i + j] >= NUMBER_OF_VECTORS_H)
                grid_check_failed();
            dualgrid -> vorticity_signs[4*i + j] = vorticity_signs[4*i + j];
            if (dualgrid -> vorticity_signs[4*i + j] != -1 && dualgrid -> vorticity_signs[4*i + j] != 1)
                grid_check_failed();
        }
    }
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        grid -> volume[i] = volume[i];
        if (grid -> volume[i] <= 0)
            grid_check_failed();
        grid -> z_scalar[i] = z_scalar[i];
        if (grid -> z_scalar[i] <= 0)
            grid_check_failed();
       grid -> gravity_potential[i] = gravity_potential[i];
       if (grid -> gravity_potential[i] <= 0)
       		grid_check_failed();
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; ++i)
    {
        grid -> normal_distance[i] = normal_distance[i];
        if (grid -> normal_distance[i] <= 0)
        	grid_check_failed();
        grid -> area[i] = area[i];
        if (grid -> area[i] <= 0)
            grid_check_failed();
        grid -> z_vector[i] = z_vector[i];
        if (grid -> z_vector[i] < -400)
            grid_check_failed();
        grid -> gravity[i] = gravity[i];
        if (grid -> gravity[i] > 0 || grid -> gravity[i] < -10)
        	grid_check_failed;
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
    {
        grid -> to_index[i] = to_index[i];
        if (grid -> to_index[i] >= NUMBER_OF_SCALARS_H || grid -> to_index[i] < 0)
            grid_check_failed();
        grid -> from_index[i] = from_index[i];
        if (grid -> from_index[i] >= NUMBER_OF_SCALARS_H || grid -> from_index[i] < 0)
            grid_check_failed();
        grid -> direction[i] = direction[i];
        if (fabs(grid -> direction[i]) >= 1.0001*M_PI)
            grid_check_failed();
        for (int j = 0; j < 2; ++j)
        {
            grid -> recov_hor_ver_dual_index[2*i + j] = recov_hor_ver_dual_index[2*i + j];
            if (grid -> recov_hor_ver_dual_index[2*i + j] >= 2*NUMBER_OF_DUAL_VECTORS_V + NUMBER_OF_DUAL_VECTORS_H || grid -> recov_hor_ver_dual_index[2*i + j] < 0)
                grid_check_failed();
            grid -> recov_hor_ver_dual_weight[2*i + j] = recov_hor_ver_dual_weight[2*i + j];
            if (fabs(grid -> recov_hor_ver_dual_weight[2*i + j]) >= 1.0001)
                grid_check_failed();
            grid -> recov_hor_par_dual_index[2*i + j] = recov_hor_par_dual_index[2*i + j];
            if (grid -> recov_hor_par_dual_index[2*i + j] >= 2*NUMBER_OF_DUAL_VECTORS_H + NUMBER_OF_DUAL_VECTORS_V || grid -> recov_hor_par_dual_index[2*i + j] < 0)
                grid_check_failed();
            grid -> recov_hor_par_dual_weight[2*i + j] = recov_hor_par_dual_weight[2*i + j];
            if (fabs(grid -> recov_hor_par_dual_weight[2*i + j]) >= 1.0001)
                grid_check_failed();
        }
        for (int j = 0; j < 10; ++j)
        {
            grid -> recov_hor_par_pri_index[10*i + j] = recov_hor_par_pri_index[10*i + j];
            if (grid -> recov_hor_par_pri_index[10*i + j] >= NUMBER_OF_VECTORS_H || grid -> recov_hor_par_pri_index[10*i + j] < 0)
                grid_check_failed();
            grid -> recov_hor_par_pri_weight[10*i + j] = recov_hor_par_pri_weight[10*i + j];
            if (fabs(grid -> recov_hor_par_pri_weight[10*i + j]) >= 0.30)
                grid_check_failed();
		}
        for (int j = 0; j < 4; ++j)
        {
            grid -> recov_hor_ver_pri_index[4*i + j] = recov_hor_ver_pri_index[4*i + j];
            if (grid -> recov_hor_ver_pri_index[4*i + j] >= 2*NUMBER_OF_VECTORS_V + NUMBER_OF_VECTORS_H || grid -> recov_hor_ver_pri_index[4*i + j] < 0)
                grid_check_failed();
        }
    }
    for (int i = 0; i < NUMBER_OF_VECTORS_V; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            grid -> recov_ver_0_pri_index[6*i + j] = recov_ver_0_pri_index[6*i + j];
            if (grid -> recov_ver_0_pri_index[6*i + j] >= 2*NUMBER_OF_VECTORS_H + NUMBER_OF_VECTORS_V || grid -> recov_ver_0_pri_index[6*i + j] < 0)
                grid_check_failed();
            grid -> recov_ver_0_pri_weight[6*i + j] = recov_ver_0_pri_weight[6*i + j];
            if (fabs(grid -> recov_ver_0_pri_weight[6*i + j]) >= 1.0001)
                grid_check_failed();
            grid -> recov_ver_0_dual_index[6*i + j] = recov_ver_0_dual_index[6*i + j];
            if (grid -> recov_ver_0_dual_index[6*i + j] >= NUMBER_OF_DUAL_VECTORS_H || grid -> recov_ver_0_dual_index[6*i + j] < 0)
                grid_check_failed();
            grid -> recov_ver_0_dual_weight[6*i + j] = recov_ver_0_dual_weight[6*i + j];
            if (fabs(grid -> recov_ver_0_dual_weight[6*i + j]) >= 1.0001)
                grid_check_failed();
            grid -> recov_ver_1_pri_index[6*i + j] = recov_ver_1_pri_index[6*i + j];
            if (recov_ver_1_pri_index[6*i + j] >= 2*NUMBER_OF_VECTORS_H + NUMBER_OF_VECTORS_V || grid -> recov_ver_1_pri_index[6*i + j] < 0)
                grid_check_failed();
            grid -> recov_ver_1_pri_weight[6*i + j] = recov_ver_1_pri_weight[6*i + j];
            if (fabs(grid -> recov_ver_1_pri_weight[6*i + j]) >= 1.0001)
                grid_check_failed();
            grid -> recov_ver_1_dual_index[6*i + j] = recov_ver_1_dual_index[6*i + j];
            if (grid -> recov_ver_1_dual_index[6*i + j] >= NUMBER_OF_DUAL_VECTORS_H || grid -> recov_ver_1_dual_index[6*i + j] < 0)
                grid_check_failed();
            grid -> recov_ver_1_dual_weight[6*i + j] = recov_ver_1_dual_weight[6*i + j];
            if (fabs(grid -> recov_ver_1_dual_weight[6*i + j]) >= 1.0001)
                grid_check_failed();
        }
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_PER_LAYER; ++i)
    {
        dualgrid -> f_vec[i] = f_vec[i];
        if (fabs(dualgrid -> f_vec[i]) > 2*OMEGA)
            grid_check_failed();
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS; ++i)
    {
        dualgrid -> normal_distance[i] = normal_distance_dual[i];
        if (dualgrid -> normal_distance[i] <= 0)
            grid_check_failed();
    }
    for (int i = 0; i < NUMBER_OF_DUAL_H_VECTORS + NUMBER_OF_H_VECTORS; ++i)
    {
        dualgrid -> area[i] = area_dual[i];
        if (dualgrid -> area[i] <= 0)
            grid_check_failed();
    }
    for (int i = 0; i < NUMBER_OF_DUAL_VECTORS_H; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            dualgrid -> h_curl_indices[4*i + j] = h_curl_indices[4*i + j];
            if (dualgrid -> h_curl_indices[4*i + j] >= 2*NUMBER_OF_VECTORS_H + NUMBER_OF_VECTORS_V || dualgrid -> h_curl_indices[4*i + j] < 0)
                grid_check_failed();
            dualgrid -> h_curl_signs[4*i + j] = h_curl_signs[4*i + j];
            if (dualgrid -> h_curl_signs[4*i + j] != -1 && dualgrid -> h_curl_signs[4*i + j] != 1)
                grid_check_failed();
        }
    }
    for (int i = 0; i < NUMBER_OF_ORO_LAYERS + 1; ++i)
    {
    	for (int j = 0; j < NUMBER_OF_VECTORS_V; ++j)
    	{
    		for (int k = 0; k < 3; ++k)
    		{
				grid -> vertical_contravar_unit[i*3*NUMBER_OF_VECTORS_V + 3*j + k] = vertical_contravar_unit[i*3*NUMBER_OF_VECTORS_V + 3*j + k];
				if (fabs(grid -> vertical_contravar_unit[i*3*NUMBER_OF_VECTORS_V + 3*j + k]) > 1.000001)
					grid_check_failed();
			}
    	}
    }
    printf("passed\n");
    free(vertical_contravar_unit);
    free(gravity);
    free(direction);
    free(normal_distance);
    free(volume);
    free(area);
    free(z_scalar);
    free(z_vector);
    free(recov_hor_par_dual_weight);
    free(recov_hor_ver_dual_weight);
    free(recov_hor_par_pri_weight);
    free(recov_ver_0_pri_weight);
    free(recov_ver_0_dual_weight);
    free(recov_ver_1_pri_weight);
    free(recov_ver_1_dual_weight);
    free(adjacent_vector_indices_dual_h);
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
    free(recov_ver_0_pri_index);
    free(recov_ver_0_dual_index);
    free(recov_ver_1_pri_index);
    free(recov_ver_1_dual_index);
    free(adjacent_signs_h);
    free(vorticity_signs);
    free(h_curl_signs);
    return 0;
}

int calc_delta_t(double cfl_margin, double *delta_t, Grid *grid)
{
    double max_speed = 350;
    double min_dist_horizontal = RADIUS;
    double min_dist_vertical = RADIUS;
    double delta_t_horizontal;
    double delta_t_vertical;
    for (int i = 0; i < NUMBER_OF_LAYERS; ++i)
    {
        for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
        {
            if (grid -> normal_distance[NUMBER_OF_VECTORS_V + i*NUMBER_OF_VECTORS_PER_LAYER + j] < min_dist_horizontal)
                min_dist_horizontal = grid -> normal_distance[NUMBER_OF_VECTORS_V + i*NUMBER_OF_VECTORS_PER_LAYER + j];
        }
    }
    for (int i = 1; i < NUMBER_OF_LEVELS - 1; ++i)
    {
        for (int j = 0; j < NUMBER_OF_VECTORS_V; ++j)
        {
            if (grid -> normal_distance[i*NUMBER_OF_VECTORS_PER_LAYER + j] < min_dist_vertical)
                min_dist_vertical = grid -> normal_distance[i*NUMBER_OF_VECTORS_PER_LAYER + j];
        }
    }
    delta_t_horizontal = (1 - cfl_margin)*min_dist_horizontal/max_speed;
    delta_t_vertical = (1 - cfl_margin)*min_dist_vertical/max_speed;
    *delta_t = delta_t_vertical;
    if (delta_t_horizontal < delta_t_vertical)
        *delta_t = delta_t_horizontal;
    return 1;
}


void grid_check_failed()
{
    printf("failed\n");
    exit(2);
}


