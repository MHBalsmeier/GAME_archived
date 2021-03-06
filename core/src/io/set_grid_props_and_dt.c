/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "../enum_and_typedefs.h"
#include "../spatial_operators/spatial_operators.h"
#include "../diagnostics/diagnostics.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
int set_grid_properties(Grid *grid, Dualgrid *dualgrid, char GEO_PROP_FILE[])
{
    double *normal_distance = malloc(NO_OF_VECTORS*sizeof(double));
    double *volume = malloc(NO_OF_SCALARS*sizeof(double));
    double *area = malloc(NO_OF_VECTORS*sizeof(double));
    double *z_scalar = malloc(NO_OF_SCALARS*sizeof(double));
    double *z_vector = malloc(NO_OF_VECTORS*sizeof(double));
    double *trsk_weights = malloc(10*NO_OF_VECTORS_H*sizeof(double));
    double *area_dual = malloc((NO_OF_DUAL_H_VECTORS + NO_OF_H_VECTORS)*sizeof(double));
    double *f_vec = malloc(2*NO_OF_VECTORS_H*sizeof(double));
    double *direction = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *gravity_potential = malloc(NO_OF_SCALARS*sizeof(double));
    double *volume_ratios = malloc(2*NO_OF_SCALARS*sizeof(double));
    double *inner_product_weights = malloc(8*NO_OF_SCALARS*sizeof(double));
    double *slope = malloc(NO_OF_VECTORS*sizeof(double));
    double *remap_horpri2hordual_vector_weights = malloc(2*NO_OF_DUAL_H_VECTORS*sizeof(double));
    double *density_to_rhombus_weights = malloc(4*NO_OF_VECTORS_H*sizeof(double));
    double *normal_distance_dual = malloc(NO_OF_DUAL_VECTORS*sizeof(double));
    double *latitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    int *from_index = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *to_index = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *from_index_dual = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *to_index_dual = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *trsk_indices = malloc(10*NO_OF_VECTORS_H*sizeof(int));
    int *trsk_modified_curl_indices = malloc(10*NO_OF_VECTORS_H*sizeof(int));
    int *adjacent_vector_indices_h = malloc(6*NO_OF_SCALARS_H*sizeof(int));
    int *adjacent_vector_indices_dual_h = malloc(3*NO_OF_DUAL_SCALARS_H*sizeof(int));
    int *vorticity_indices = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *adjacent_signs_h = malloc(6*NO_OF_SCALARS_H*sizeof(int));
    int *vorticity_signs = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *density_to_rhombus_indices = malloc(4*NO_OF_VECTORS_H*sizeof(int));
    int *no_of_shaded_points_scalar = malloc(NO_OF_SCALARS_H*sizeof(int));
    int *no_of_shaded_points_vector = malloc(NO_OF_VECTORS_H*sizeof(int));
    int ncid, retval;
    int normal_distance_id, volume_id, area_id, z_scalar_id, z_vector_id, trsk_weights_id, area_dual_id, f_vec_id, to_index_id, from_index_id, to_index_dual_id, from_index_dual_id, adjacent_vector_indices_h_id, vorticity_indices_id, trsk_indices_id, trsk_modified_curl_indices_id, adjacent_signs_h_id, vorticity_signs_id, direction_id, gravity_potential_id, inner_product_weights_id, slope_id, volume_ratios_id, remap_horpri2hordual_vector_weights_id, density_to_rhombus_weights_id, density_to_rhombus_indices_id, normal_distance_dual_id, adjacent_vector_indices_dual_h_id, latitude_scalar_id, longitude_scalar_id, stretching_parameter_id, no_of_shaded_points_scalar_id, no_of_shaded_points_vector_id;
    double stretching_parameter;
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
    if ((retval = nc_inq_varid(ncid, "trsk_weights", &trsk_weights_id)))
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
    if ((retval = nc_inq_varid(ncid, "trsk_indices", &trsk_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "trsk_modified_curl_indices", &trsk_modified_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_signs_h", &adjacent_signs_h_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vorticity_signs", &vorticity_signs_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "inner_product_weights", &inner_product_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "slope", &slope_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "remap_horpri2hordual_vector_weights", &remap_horpri2hordual_vector_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_to_rhombus_weights", &density_to_rhombus_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_to_rhombus_indices", &density_to_rhombus_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "latitude_scalar", &latitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_scalar", &longitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "stretching_parameter", &stretching_parameter_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "no_of_shaded_points_scalar", &no_of_shaded_points_scalar_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "no_of_shaded_points_vector", &no_of_shaded_points_vector_id)))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, stretching_parameter_id, &stretching_parameter)))
        ERR(retval);
    grid -> stretching_parameter = stretching_parameter;
    if ((retval = nc_get_var_double(ncid, normal_distance_id, &normal_distance[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, inner_product_weights_id, &inner_product_weights[0])))
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
    if ((retval = nc_get_var_double(ncid, trsk_weights_id, &trsk_weights[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, area_dual_id, &area_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, direction_id, &direction[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, slope_id, &slope[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, f_vec_id, &f_vec[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, remap_horpri2hordual_vector_weights_id, &remap_horpri2hordual_vector_weights[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, density_to_rhombus_weights_id, &density_to_rhombus_weights[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, normal_distance_dual_id, &normal_distance_dual[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_scalar_id, &latitude_scalar[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_scalar_id, &longitude_scalar[0])))
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
    if ((retval = nc_get_var_int(ncid, trsk_indices_id, &trsk_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, trsk_modified_curl_indices_id, &trsk_modified_curl_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, adjacent_signs_h_id, &adjacent_signs_h[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, vorticity_signs_id, &vorticity_signs[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, density_to_rhombus_indices_id, &density_to_rhombus_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, no_of_shaded_points_scalar_id, &no_of_shaded_points_scalar[0])))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, no_of_shaded_points_vector_id, &no_of_shaded_points_vector[0])))
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
            }
            grid -> adjacent_signs_h[6*i + j] = adjacent_signs_h[6*i + j];
        }
        grid -> latitude_scalar[i] = latitude_scalar[i];
        grid -> longitude_scalar[i] = longitude_scalar[i];
        grid -> no_of_shaded_points_scalar[i] = no_of_shaded_points_scalar[i];
    }
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            dualgrid -> vorticity_indices[4*i + j] = vorticity_indices[4*i +j];
            dualgrid -> vorticity_signs[4*i + j] = vorticity_signs[4*i + j];
            grid -> density_to_rhombus_indices[4*i + j] = density_to_rhombus_indices[4*i + j];
            grid -> density_to_rhombus_weights[4*i + j] = density_to_rhombus_weights[4*i + j];
        }
        grid -> to_index[i] = to_index[i];
        grid -> from_index[i] = from_index[i];
        grid -> direction[i] = direction[i];
        dualgrid -> from_index[i] = from_index[i];
        dualgrid -> to_index[i] = to_index[i];
        grid -> no_of_shaded_points_vector[i] = no_of_shaded_points_vector[i];
        for (int j = 0; j < 10; ++j)
        {
            grid -> trsk_indices[10*i + j] = trsk_indices[10*i + j];
            grid -> trsk_modified_curl_indices[10*i + j] = trsk_modified_curl_indices[10*i + j];
            grid -> trsk_weights[10*i + j] = trsk_weights[10*i + j];
		}
        for (int j = 0; j < 2; ++j)
        {
		    dualgrid -> f_vec[2*i + j] = f_vec[2*i + j];
        }
    }
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        grid -> volume[i] = volume[i];
        grid -> z_scalar[i] = z_scalar[i];
       	grid -> gravity_potential[i] = gravity_potential[i];
       	for (int j = 0; j < 8; ++j)
		{
           grid -> inner_product_weights[8*i + j] = inner_product_weights[8*i + j];
       	}
   	   	for (int j = 0; j < 2; ++j)
   	   	{
   	   		grid -> volume_ratios[2*i + j] = volume_ratios[2*i + j];
   	   	}
    }
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        grid -> normal_distance[i] = normal_distance[i];
        grid -> area[i] = area[i];
        grid -> z_vector[i] = z_vector[i];
        grid -> slope[i] = slope[i];
    }
    for (int i = 0; i < NO_OF_DUAL_VECTORS; ++i)
    {
        dualgrid -> normal_distance[i] = normal_distance_dual[i];
    }
    for (int i = 0; i < NO_OF_DUAL_H_VECTORS + NO_OF_H_VECTORS; ++i)
    {
        dualgrid -> area[i] = area_dual[i];
    }
    for (int i = 0; i < NO_OF_DUAL_H_VECTORS; ++i)
    {
    	for (int j = 0; j < 2; ++j)
    	{
    		grid -> remap_horpri2hordual_vector_weights[2*i + j] = remap_horpri2hordual_vector_weights[2*i + j];
    	}
    }
    for (int i = 0; i < NO_OF_DUAL_SCALARS_H; ++i)
    {
    	for (int j = 0; j < 3; ++j)
    	{
		    dualgrid -> adjacent_vector_indices_h[3*i + j] = adjacent_vector_indices_dual_h[3*i + j];
        }
    }
    grad(grid -> gravity_potential, grid -> gravity_m, grid);
    printf("stretching parameter of the vertical grid: %lf\n", stretching_parameter);
    free(no_of_shaded_points_scalar);
    free(no_of_shaded_points_vector);
    free(latitude_scalar);
    free(longitude_scalar);
    free(normal_distance_dual);
    free(density_to_rhombus_indices);
    free(density_to_rhombus_weights);
    free(remap_horpri2hordual_vector_weights);
    free(volume_ratios);
    free(slope);
    free(inner_product_weights);
    free(gravity_potential);
    free(direction);
    free(normal_distance);
    free(volume);
    free(area);
    free(z_scalar);
    free(z_vector);
    free(trsk_weights);
    free(area_dual);
    free(f_vec);
    free(from_index);
    free(to_index);
    free(from_index_dual);
    free(to_index_dual);
    free(adjacent_vector_indices_h);
    free(vorticity_indices);
    free(trsk_indices);
    free(trsk_modified_curl_indices);
    free(adjacent_signs_h);
    free(vorticity_signs);
    free(adjacent_vector_indices_dual_h);
    return 0;
}

int calc_delta_t_and_related(double cfl_margin, double *delta_t, Grid *grid, Dualgrid *dualgrid, State *state, Config_info *config_info)
{
    double max_sound_speed = 0;
    double sound_speed_value;
    
    // the gravity wave criterion
    double n_freq = 0.03;
    double delta_t_brunt_vaisala = 2/n_freq;
    // delta_t_brunt_vaisala = (1 - cfl_margin)*delta_t_brunt_vaisala;
    
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	sound_speed_value = pow(gas_constant_diagnostics(state, i, config_info)
    	*spec_heat_cap_diagnostics_p(state, i, config_info)/spec_heat_cap_diagnostics_v(state, i, config_info)
    	*state -> temperature_gas[i], 0.5);
    	if (sound_speed_value > max_sound_speed)
    	{
    		max_sound_speed = sound_speed_value;
    	}
    }
    // adding a safety margin
    max_sound_speed = 1.1*max_sound_speed;
    double min_dist_horizontal = RADIUS;
    double edge_area = 0;
    for (int i = 0; i < NO_OF_LAYERS; ++i)
    {
        for (int j = 0; j < NO_OF_VECTORS_H; ++j)
        {
            edge_area += 0.5*grid -> normal_distance[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j]
            *dualgrid -> normal_distance[i*NO_OF_DUAL_VECTORS_PER_LAYER + j];
            if (grid -> normal_distance[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] < min_dist_horizontal)
            {
                min_dist_horizontal = grid -> normal_distance[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j];
    		}
        }
    }
	*delta_t = (1 - cfl_margin)*min_dist_horizontal/max_sound_speed;
	if (*delta_t > delta_t_brunt_vaisala)
	{
		*delta_t = delta_t_brunt_vaisala;
	}
	
	grid -> mean_area_edge = edge_area/NO_OF_H_VECTORS;
	
	/*
	the homogeneous prefactor of the divergence damping
	this is for RES_ID = 5
	double div_damp_coeff = 0.5e15; // unstable
	div_damp_coeff = 5e15; // grid-scale noise is produced, leading to instability
	div_damp_coeff = 28e15; // grid-scale noise is produced, leading to instability (best choice probably)
	div_damp_coeff = 250e15; // too close to instability, unrealistically large compared to the literature
	div_damp_coeff = 2500e15; // unstable
	*/
	// generalized version following the ICON paper
	config_info -> div_damp_coeff = 1/(500*(*delta_t))*pow(grid -> mean_area_edge, 2);
    return 1;
}



