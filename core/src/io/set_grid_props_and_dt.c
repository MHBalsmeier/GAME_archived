/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
This file contains functions for reading the grid properties as well as setting the time step.
*/

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "geos95.h"
#include "../enum_and_typedefs.h"
#include "../spatial_operators/spatial_operators.h"
#include "../thermodynamics/thermodynamics.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int set_grid_properties(Grid *grid, Dualgrid *dualgrid, char GEO_PROP_FILE[])
{
    int ncid, retval;
    int normal_distance_id, volume_id, area_id, z_scalar_id, z_vector_id, trsk_weights_id, area_dual_id, z_vector_dual_id, f_vec_id, to_index_id, from_index_id, to_index_dual_id, from_index_dual_id, adjacent_vector_indices_h_id, trsk_indices_id, trsk_modified_curl_indices_id, adjacent_signs_h_id, direction_id, gravity_potential_id, inner_product_weights_id, density_to_rhombi_weights_id, density_to_rhombi_indices_id, normal_distance_dual_id, vorticity_indices_triangles_id, vorticity_signs_triangles_id, latitude_scalar_id, longitude_scalar_id, stretching_parameter_id, no_of_shaded_points_scalar_id, no_of_shaded_points_vector_id, interpol_indices_id, interpol_weights_id;
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
    if ((retval = nc_inq_varid(ncid, "trsk_weights", &trsk_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "area_dual", &area_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_vector_dual", &z_vector_dual_id)))
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
    if ((retval = nc_inq_varid(ncid, "vorticity_indices_triangles", &vorticity_indices_triangles_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vorticity_signs_triangles", &vorticity_signs_triangles_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "trsk_indices", &trsk_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "trsk_modified_curl_indices", &trsk_modified_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_signs_h", &adjacent_signs_h_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "inner_product_weights", &inner_product_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_to_rhombi_weights", &density_to_rhombi_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_to_rhombi_indices", &density_to_rhombi_indices_id)))
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
    if ((retval = nc_inq_varid(ncid, "interpol_indices", &interpol_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "interpol_weights", &interpol_weights_id)))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, stretching_parameter_id, &stretching_parameter)))
        ERR(retval);
    grid -> stretching_parameter = stretching_parameter;
    if ((retval = nc_get_var_double(ncid, normal_distance_id, &(grid -> normal_distance[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, inner_product_weights_id, &(grid -> inner_product_weights[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, volume_id, &(grid -> volume[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, area_id, &(grid -> area[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_scalar_id, &(grid -> z_scalar[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, gravity_potential_id, &(grid -> gravity_potential[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_vector_id, &(grid -> z_vector[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, trsk_weights_id, &(grid -> trsk_weights[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, area_dual_id, &(dualgrid -> area[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_vector_dual_id, &(dualgrid -> z_vector[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, direction_id, &(grid -> direction[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, f_vec_id, &(dualgrid -> f_vec[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, density_to_rhombi_weights_id, &(grid -> density_to_rhombi_weights[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, normal_distance_dual_id, &(dualgrid -> normal_distance[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_scalar_id, &(grid -> latitude_scalar[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_scalar_id, &(grid -> longitude_scalar[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, interpol_weights_id, &(grid -> latlon_interpol_weights[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, from_index_id, &(grid -> from_index[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, to_index_id, &(grid -> to_index[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, from_index_dual_id, &(dualgrid -> from_index[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, to_index_dual_id, &(dualgrid -> to_index[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, adjacent_vector_indices_h_id, &(grid -> adjacent_vector_indices_h[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, vorticity_indices_triangles_id, &(dualgrid -> vorticity_indices_triangles[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, vorticity_signs_triangles_id, &(dualgrid -> vorticity_signs_triangles[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, trsk_indices_id, &(grid -> trsk_indices[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, trsk_modified_curl_indices_id, &(grid -> trsk_modified_curl_indices[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, adjacent_signs_h_id, &(grid -> adjacent_signs_h[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, density_to_rhombi_indices_id, &(grid -> density_to_rhombi_indices[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, no_of_shaded_points_scalar_id, &(grid -> no_of_shaded_points_scalar[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, no_of_shaded_points_vector_id, &(grid -> no_of_shaded_points_vector[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, interpol_indices_id, &(grid -> latlon_interpol_indices[0]))))
        ERR(retval);
    if ((retval = nc_close(ncid)))
        ERR(retval);
    for (int i = 0; i < 6*NO_OF_SCALARS_H; ++i)
    {
        if (grid -> adjacent_vector_indices_h[i] == -1)
        {
        	grid -> adjacent_vector_indices_h[i] = 0;
        }
    }
    // determining coordinate slopes
    grad_hor_cov(grid -> z_scalar, grid -> slope, grid);
    // computing the gradient of the gravity potential
    grad(grid -> gravity_potential, grid -> gravity_m, grid);
    printf("stretching parameter of the vertical grid: %lf\n", stretching_parameter);
    return 0;
}

int calc_delta_t_and_related(double cfl_margin, double *delta_t, Grid *grid, Dualgrid *dualgrid, State *state, Config_info *config_info)
{
	/*
	This function sets the timestep of the model.
	*/
	
    double max_sound_speed = 0;
    double sound_speed_value;
    
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
    for (int i = 0; i < NO_OF_LAYERS; ++i)
    {
        for (int j = 0; j < NO_OF_VECTORS_H; ++j)
        {
            if (grid -> normal_distance[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] < min_dist_horizontal)
            {
                min_dist_horizontal = grid -> normal_distance[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j];
    		}
        }
    }
	*delta_t = (1 - cfl_margin)*min_dist_horizontal/max_sound_speed;
	
	// diffusion preparation
	int layer_index, h_index;
	double cell_area_sum = 0;
	for (int i = 0; i < NO_OF_LEVELS*NO_OF_SCALARS_H; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		cell_area_sum += grid -> area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
	}
	grid -> mean_area_cell = cell_area_sum/(NO_OF_LEVELS*NO_OF_SCALARS_H);
    return 1;
}










