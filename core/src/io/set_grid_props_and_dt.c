/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "../enum_and_typedefs.h"
#include "../spatial_operators/spatial_operators.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

void grid_check_failed();

int set_grid_properties(Grid *grid, Dualgrid *dualgrid, char GEO_PROP_FILE[])
{
    double *normal_distance = malloc(NO_OF_VECTORS*sizeof(double));
    double *volume = malloc(NO_OF_SCALARS*sizeof(double));
    double *area = malloc(NO_OF_VECTORS*sizeof(double));
    double *z_scalar = malloc(NO_OF_SCALARS*sizeof(double));
    double *z_vector = malloc(NO_OF_VECTORS*sizeof(double));
    double *trsk_modified_weights = malloc(10*NO_OF_VECTORS_H*sizeof(double));
    double *recov_ver_weight = malloc(6*NO_OF_LEVELS*NO_OF_SCALARS_H*sizeof(double));
    double *area_dual = malloc((NO_OF_DUAL_H_VECTORS + NO_OF_H_VECTORS)*sizeof(double));
    double *f_vec = malloc(3*NO_OF_VECTORS_H*sizeof(double));
    double *direction = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *gravity_potential = malloc(NO_OF_SCALARS*sizeof(double));
    double *volume_ratios = malloc(2*NO_OF_SCALARS*sizeof(double));
    double *e_kin_weights = malloc(8*NO_OF_SCALARS*sizeof(double));
    double *slope = malloc(NO_OF_VECTORS*sizeof(double));
    double *recov_primal2dual_weights = malloc(2*NO_OF_DUAL_H_VECTORS*sizeof(double));
    double *density_to_rhombus_weights = malloc(4*NO_OF_VECTORS_H*sizeof(double));
    double *normal_distance_dual = malloc(NO_OF_DUAL_VECTORS*sizeof(double));
    int *from_index = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *to_index = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *from_index_dual = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *to_index_dual = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *trsk_modified_velocity_indices = malloc(10*NO_OF_VECTORS_H*sizeof(int));
    int *trsk_modified_curl_indices = malloc(10*NO_OF_VECTORS_H*sizeof(int));
    int *adjacent_vector_indices_h = malloc(6*NO_OF_SCALARS_H*sizeof(int));
    int *adjacent_vector_indices_dual_h = malloc(3*NO_OF_DUAL_SCALARS_H*sizeof(int));
    int *vorticity_indices = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *h_curl_indices = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *adjacent_signs_h = malloc(6*NO_OF_SCALARS_H*sizeof(int));
    int *vorticity_signs = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *h_curl_signs = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *density_to_rhombus_indices = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int ncid, retval;
    int normal_distance_id, volume_id, area_id, z_scalar_id, z_vector_id, trsk_modified_weights_id, recov_ver_weight_id, area_dual_id, f_vec_id, to_index_id, from_index_id, to_index_dual_id, from_index_dual_id, adjacent_vector_indices_h_id, vorticity_indices_id, h_curl_indices_id, trsk_modified_velocity_indices_id, trsk_modified_curl_indices_id, adjacent_signs_h_id, vorticity_signs_id, h_curl_signs_id, direction_id, gravity_potential_id, e_kin_weights_id, slope_id, volume_ratios_id, recov_primal2dual_weights_id, density_to_rhombus_weights_id, density_to_rhombus_indices_id, normal_distance_dual_id, adjacent_vector_indices_dual_h_id;
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
    if ((retval = nc_inq_varid(ncid, "z_vector", &z_vector_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "volume_ratios", &volume_ratios_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "trsk_modified_weights", &trsk_modified_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_ver_weight", &recov_ver_weight_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "area_dual", &area_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "f_vec", &f_vec_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "to_index", &to_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "to_index_dual", &to_index_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "direction", &direction_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "normal_distance_dual", &normal_distance_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "from_index", &from_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "from_index_dual", &from_index_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_vector_indices_h", &adjacent_vector_indices_h_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_vector_indices_dual_h", &adjacent_vector_indices_dual_h_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vorticity_indices", &vorticity_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "h_curl_indices", &h_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "trsk_modified_velocity_indices", &trsk_modified_velocity_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "trsk_modified_curl_indices", &trsk_modified_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_signs_h", &adjacent_signs_h_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vorticity_signs", &vorticity_signs_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "h_curl_signs", &h_curl_signs_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "e_kin_weights", &e_kin_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "slope", &slope_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "recov_primal2dual_weights", &recov_primal2dual_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_to_rhombus_weights", &density_to_rhombus_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_to_rhombus_indices", &density_to_rhombus_indices_id)))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, normal_distance_id, &normal_distance[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, e_kin_weights_id, &e_kin_weights[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, volume_id, &volume[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, area_id, &area[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_scalar_id, &z_scalar[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, gravity_potential_id, &gravity_potential[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_vector_id, &z_vector[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, volume_ratios_id, &volume_ratios[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, trsk_modified_weights_id, &trsk_modified_weights[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_ver_weight_id, &recov_ver_weight[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, area_dual_id, &area_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, direction_id, &direction[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, slope_id, &slope[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, f_vec_id, &f_vec[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, recov_primal2dual_weights_id, &recov_primal2dual_weights[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, density_to_rhombus_weights_id, &density_to_rhombus_weights[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, normal_distance_dual_id, &normal_distance_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, from_index_id, &from_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, to_index_id, &to_index[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, adjacent_vector_indices_h_id, &adjacent_vector_indices_h[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, adjacent_vector_indices_dual_h_id, &adjacent_vector_indices_dual_h[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, vorticity_indices_id, &vorticity_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, h_curl_indices_id, &h_curl_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, trsk_modified_velocity_indices_id, &trsk_modified_velocity_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, trsk_modified_curl_indices_id, &trsk_modified_curl_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, adjacent_signs_h_id, &adjacent_signs_h[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, vorticity_signs_id, &vorticity_signs[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, h_curl_signs_id, &h_curl_signs[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, density_to_rhombus_indices_id, &density_to_rhombus_indices[0])))
        ERR(retval);
    if ((retval = nc_close(ncid)))
        ERR(retval);
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            grid -> adjacent_vector_indices_h[6*i + j] = adjacent_vector_indices_h[6*i + j];
            if (grid -> adjacent_vector_indices_h[6*i + j] >= NO_OF_VECTORS_H || grid -> adjacent_vector_indices_h[6*i + j] < 0)
            {
                if (grid -> adjacent_vector_indices_h[6*i + j] == -1)
                {
                	grid -> adjacent_vector_indices_h[6*i + j] = 0;
            	}
            	else
            	{
                	grid_check_failed();
            	}
            }
            grid -> adjacent_signs_h[6*i + j] = adjacent_signs_h[6*i + j];
            if (grid -> adjacent_signs_h[6*i + j] != -1 && grid -> adjacent_signs_h[6*i + j] != 1)
            {
            	if (i >= NO_OF_PENTAGONS || j != 5 || grid -> adjacent_signs_h[6*i + j] != 0)
            	{
               		grid_check_failed(); 	
            	}
        	}
        }
    }
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            dualgrid -> vorticity_indices[4*i + j] = vorticity_indices[4*i +j];
            if (dualgrid -> vorticity_indices[4*i + j] >= NO_OF_VECTORS_H)
                grid_check_failed();
            dualgrid -> vorticity_signs[4*i + j] = vorticity_signs[4*i + j];
            if (dualgrid -> vorticity_signs[4*i + j] != -1 && dualgrid -> vorticity_signs[4*i + j] != 1)
                grid_check_failed();
            dualgrid -> h_curl_indices[4*i + j] = h_curl_indices[4*i + j];
            if (dualgrid -> h_curl_indices[4*i + j] >= 2*NO_OF_VECTORS_H + NO_OF_SCALARS_H || dualgrid -> h_curl_indices[4*i + j] < 0)
                grid_check_failed();
            dualgrid -> h_curl_signs[4*i + j] = h_curl_signs[4*i + j];
            if (dualgrid -> h_curl_signs[4*i + j] != -1 && dualgrid -> h_curl_signs[4*i + j] != 1)
                grid_check_failed();
            grid -> density_to_rhombus_indices[4*i + j] = density_to_rhombus_indices[4*i + j];
            if (grid -> density_to_rhombus_indices[4*i + j] < 0 && grid -> density_to_rhombus_indices[4*i + j] >= NO_OF_VECTORS_H)
                grid_check_failed();
            grid -> density_to_rhombus_weights[4*i + j] = density_to_rhombus_weights[4*i + j];
            if (grid -> density_to_rhombus_weights[4*i + j] < 0 && grid -> density_to_rhombus_weights[4*i + j] >= 1)
                grid_check_failed();
        }
        grid -> to_index[i] = to_index[i];
        if (grid -> to_index[i] >= NO_OF_SCALARS_H || grid -> to_index[i] < 0)
            grid_check_failed();
        grid -> from_index[i] = from_index[i];
        if (grid -> from_index[i] >= NO_OF_SCALARS_H || grid -> from_index[i] < 0)
            grid_check_failed();
        grid -> direction[i] = direction[i];
        if (fabs(grid -> direction[i]) >= 1.0001*M_PI)
            grid_check_failed();
        dualgrid -> from_index[i] = from_index[i];
        if (dualgrid -> from_index[i] < 0 || dualgrid -> from_index[i] >= NO_OF_DUAL_SCALARS_H)
        	grid_check_failed();
        dualgrid -> to_index[i] = to_index[i];
        if (dualgrid -> to_index[i] < 0 || dualgrid -> to_index[i] >= NO_OF_DUAL_SCALARS_H)
        	grid_check_failed();
        for (int j = 0; j < 10; ++j)
        {
            grid -> trsk_modified_velocity_indices[10*i + j] = trsk_modified_velocity_indices[10*i + j];
            if (grid -> trsk_modified_velocity_indices[10*i + j] >= NO_OF_VECTORS_H || grid -> trsk_modified_velocity_indices[10*i + j] < 0)
                grid_check_failed();
            grid -> trsk_modified_curl_indices[10*i + j] = trsk_modified_curl_indices[10*i + j];
            if (grid -> trsk_modified_curl_indices[10*i + j] >= NO_OF_VECTORS_H|| grid -> trsk_modified_curl_indices[10*i + j] < 0)
                grid_check_failed();
            grid -> trsk_modified_weights[10*i + j] = trsk_modified_weights[10*i + j];
            if (fabs(grid -> trsk_modified_weights[10*i + j]) >= 0.30)
                grid_check_failed();
		}
        for (int j = 0; j < 3; ++j)
        {
		    dualgrid -> f_vec[3*i + j] = f_vec[3*i + j];
		    if (fabs(dualgrid -> f_vec[3*i + j]) > 2*OMEGA)
		        grid_check_failed();
        }
    }
    double check_sum = 0;
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        grid -> volume[i] = volume[i];
        if (grid -> volume[i] <= 0)
        {
            grid_check_failed();   
      	}
        grid -> z_scalar[i] = z_scalar[i];
        if (grid -> z_scalar[i] <= 0)
        {
            grid_check_failed();   
      	}
       	grid -> gravity_potential[i] = gravity_potential[i];
       	if (grid -> gravity_potential[i] <= 0)
        {
            grid_check_failed();   
      	}
        if (grid -> z_scalar[i] <= 0)
        {
            grid_check_failed();   
  		}
       	check_sum = 0;
       	for (int j = 0; j < 8; ++j)
		{
           grid -> e_kin_weights[8*i + j] = e_kin_weights[8*i + j];
           if (grid -> e_kin_weights[8*i + j] > 0.3 || grid -> e_kin_weights[8*i + j] < 0)
           	   grid_check_failed();
           check_sum += grid -> e_kin_weights[8*i + j];
       	}
       	if (fabs(check_sum - 1.5) > 1e-10)
       	{
			grid_check_failed();
   	   	}
   	   	check_sum = 0;
   	   	for (int j = 0; j < 2; ++j)
   	   	{
   	   		grid -> volume_ratios[2*i + j] = volume_ratios[2*i + j];
   			check_sum += grid -> volume_ratios[2*i + j];
   	   	}
   		if (fabs(check_sum - 1) > 1e-10)
   		{
   			grid_check_failed();
   		}
    }
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        grid -> normal_distance[i] = normal_distance[i];
        if (grid -> normal_distance[i] <= 0)
        	grid_check_failed();
        grid -> area[i] = area[i];
        if (grid -> area[i] <= 0)
            grid_check_failed();
        grid -> z_vector[i] = z_vector[i];
        if (grid -> z_vector[i] < -600)
            grid_check_failed();
        grid -> slope[i] = slope[i];
        if (fabs(grid -> slope[i]) > 1)
        	grid_check_failed();
    }
    for (int i = 0; i < NO_OF_DUAL_VECTORS; ++i)
    {
        dualgrid -> normal_distance[i] = normal_distance_dual[i];
        if (dualgrid -> normal_distance[i] <= 0)
        	grid_check_failed();
    }
    for (int i = 0; i < NO_OF_LEVELS*NO_OF_SCALARS_H; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            grid -> recov_ver_weight[6*i + j] = recov_ver_weight[6*i + j];
            if (fabs(grid -> recov_ver_weight[6*i + j]) >= 1.0001)
                grid_check_failed();
        }
    }
    for (int i = 0; i < NO_OF_DUAL_H_VECTORS + NO_OF_H_VECTORS; ++i)
    {
        dualgrid -> area[i] = area_dual[i];
        if (dualgrid -> area[i] <= 0)
            grid_check_failed();
    }
    for (int i = 0; i < NO_OF_DUAL_H_VECTORS; ++i)
    {
    	check_sum = 0;
    	for (int j = 0; j < 2; ++j)
    	{
    		grid -> recov_primal2dual_weights[2*i + j] = recov_primal2dual_weights[2*i + j];
    		check_sum += grid -> recov_primal2dual_weights[2*i + j];
    		if (grid -> recov_primal2dual_weights[2*i + j] < 0 || grid -> recov_primal2dual_weights[2*i + j] >= 1)
    		{
    			grid_check_failed();
    		}
    	}
		if (fabs(check_sum - 1) > 1e-10 && i >= NO_OF_VECTORS_H && i < NO_OF_DUAL_H_VECTORS - NO_OF_VECTORS_H)
		{
			grid_check_failed();
		}
    }
    for (int i = 0; i < NO_OF_DUAL_SCALARS_H; ++i)
    {
    	for (int j = 0; j < 3; ++j)
    	{
		    dualgrid -> adjacent_vector_indices_h[3*i + j] = adjacent_vector_indices_dual_h[3*i + j];
		    if (dualgrid -> adjacent_vector_indices_h[3*i + j] < 0 || dualgrid -> adjacent_vector_indices_h[3*i + j] >= NO_OF_VECTORS_H)
		        grid_check_failed();
        }
    }
    grad(grid -> gravity_potential, grid -> gravity_m, grid);
    printf("passed\n");
    free(normal_distance_dual);
    free(density_to_rhombus_indices);
    free(density_to_rhombus_weights);
    free(recov_primal2dual_weights);
    free(volume_ratios);
    free(slope);
    free(e_kin_weights);
    free(gravity_potential);
    free(direction);
    free(normal_distance);
    free(volume);
    free(area);
    free(z_scalar);
    free(z_vector);
    free(trsk_modified_weights);
    free(recov_ver_weight);
    free(area_dual);
    free(f_vec);
    free(from_index);
    free(to_index);
    free(from_index_dual);
    free(to_index_dual);
    free(adjacent_vector_indices_h);
    free(vorticity_indices);
    free(h_curl_indices);
    free(trsk_modified_velocity_indices);
    free(trsk_modified_curl_indices);
    free(adjacent_signs_h);
    free(vorticity_signs);
    free(h_curl_signs);
    free(adjacent_vector_indices_dual_h);
    return 0;
}

int calc_delta_t(double cfl_margin, double *delta_t, Grid *grid)
{
    double max_speed = 350;
    double min_dist_horizontal = RADIUS;
    double delta_t_candidate;
    double brunt_vaisala_max = 0.03;
    double delta_t_buoyancy = 2/brunt_vaisala_max;
    for (int i = 0; i < NO_OF_LAYERS; ++i)
    {
        for (int j = 0; j < NO_OF_VECTORS_H; ++j)
        {
            if (grid -> normal_distance[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] < min_dist_horizontal)
                min_dist_horizontal = grid -> normal_distance[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j];
        }
    }
    delta_t_candidate = (1 - cfl_margin)*min_dist_horizontal/max_speed;
    if (delta_t_buoyancy < delta_t_candidate)
        *delta_t = (1 - cfl_margin)*delta_t_buoyancy;
    else
    {
    	*delta_t = (1 - cfl_margin)*delta_t_candidate;
	}
    return 1;
}


void grid_check_failed()
{
    printf("failed\n");
    exit(0);
}


