/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

/*
Here, the output is written to grib and/or netcdf files and integrals are written to text files if configured that way.
*/

#include "../../../core/src/enum_and_typedefs.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <netcdf.h>
#include "io.h"
#include "../diagnostics/diagnostics.h"
#include "../spatial_operators/spatial_operators.h"
#include "../settings.h"
#include "eccodes.h"
#include "geos95.h"
#include "atmostracers.h"
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

double calc_std_dev(double [], int);

int write_out(State *state_write_out, double wind_h_lowest_layer_array[], int min_no_of_output_steps, double t_init, double t_write, Diagnostics *diagnostics, Forcings *forcings, Grid *grid, Dualgrid *dualgrid, char RUN_ID[], Io_config *io_config)
{
	// Diagnostics and forcings are primarily handed over for checks.
	
	// Time stuff.
	int write_out_divv_h;
	ask_for_divergence_output(&write_out_divv_h);
    time_t t_init_tt = (time_t) t_init;
    struct tm *p_init_time = gmtime(&t_init_tt);
    int init_year = p_init_time -> tm_year;
    int init_month = p_init_time -> tm_mon;
    int init_day = p_init_time -> tm_mday;
    int init_hour = p_init_time -> tm_hour;
    long data_date = 10000*init_year + 100*init_month + init_day;
    long data_time = init_hour;
    
    // Diagnostics of surface speicific quantities. These are always written out.
    Scalar_field *pot_temperature = calloc(1, sizeof(Scalar_field));
	double *mslp = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *surface_p = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *t2 = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *tcdc = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *rprate = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *sprate = malloc(NO_OF_SCALARS_H*sizeof(double));
	double *cape = malloc(NO_OF_SCALARS_H*sizeof(double));
    double temp_lowest_layer, pressure_value, mslp_factor, surface_p_factor, temp_mslp, temp_surface, z_height, theta_v_prime, theta_v, cape_integrand, delta_z, temp_upper, temp_lower, delta_z_temp, temperature_gradient, density_v, density_h;
	double z_tropopause = 15e3;
    double standard_vert_lapse_rate = 0.0065;
    int layer_index;
	pot_temp_diagnostics(state_write_out, *pot_temperature);
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
    	// Now the aim is to determine the value of the MSLP.
        temp_lowest_layer = state_write_out -> temperature_gas[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i];
        pressure_value = state_write_out -> density_dry[i + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H]*R_D*temp_lowest_layer + state_write_out -> tracer_densities[2*NO_OF_SCALARS + i + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H]*R_V*temp_lowest_layer;
        temp_mslp = temp_lowest_layer + standard_vert_lapse_rate*grid -> z_scalar[i + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H];
        mslp_factor = pow(1 - (temp_mslp - temp_lowest_layer)/temp_mslp, grid -> gravity_m[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i]/(R_D*standard_vert_lapse_rate));
        mslp[i] = pressure_value/mslp_factor;
        
		// Now the aim is to determine the value of the surface pressure.
		temp_surface = temp_lowest_layer + standard_vert_lapse_rate*(grid -> z_scalar[i + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + i]);
        surface_p_factor = pow(1 - (temp_surface - temp_lowest_layer)/temp_surface, grid -> gravity_m[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i]/(R_D*standard_vert_lapse_rate));
		surface_p[i] = pressure_value/surface_p_factor;
		
		// Now the aim is to calculate the 2 m temperature.
        delta_z_temp = 2 - grid -> z_scalar[i + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H];
        temp_upper = state_write_out -> temperature_gas[(NO_OF_LAYERS - 2)*NO_OF_SCALARS_H + i];
        temp_lower = temp_lowest_layer;
        temperature_gradient = (temp_upper - temp_lower)/(grid -> z_scalar[i + (NO_OF_LAYERS - 2)*NO_OF_SCALARS_H] - grid -> z_scalar[i + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H]);
        
        // Finally the temperature in 2 m height AGL can bo obtained via linear extrapolation.
        t2[i] = temp_lowest_layer + delta_z_temp*temperature_gradient;
        z_height = grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i];
        cape[i] = 0;
        density_v = state_write_out -> tracer_densities[2*NO_OF_SCALARS + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i];
        density_h = state_write_out -> density_dry[i + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H] + state_write_out -> tracer_densities[2*NO_OF_SCALARS + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i];
        theta_v_prime = (*pot_temperature)[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i]*(1 + density_v/density_h*(1/EPSILON - 1));
        layer_index = NO_OF_LAYERS - 1;
        cape[i] = 0;
        while (z_height < z_tropopause)
        {
		    density_v = state_write_out -> tracer_densities[2*NO_OF_SCALARS + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i];
		    density_h = state_write_out -> density_dry[i + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H] + state_write_out -> tracer_densities[2*NO_OF_SCALARS + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i];
            theta_v = (*pot_temperature)[layer_index*NO_OF_SCALARS_H + i]*(1 + density_v/density_h*(1/EPSILON - 1));
        	delta_z = grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + i];
        	z_height += delta_z;
        	cape_integrand = grid -> gravity_m[(NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + i]*(theta_v_prime - theta_v)/theta_v;
        	if (cape_integrand > 0)
        		cape[i] += cape_integrand*delta_z;
        	--layer_index;
        }
        
        // Now come the hydrometeors.
        sprate[i] = 0;
        rprate[i] = 0;
        tcdc[i] = 0;
        for (int k = 0; k < NO_OF_CONDENSATED_TRACERS; ++k)
        {
            for (int l = 0; l < NO_OF_LAYERS; ++l)
            {
                if (state_write_out -> tracer_densities[k*NO_OF_SCALARS + l*NO_OF_SCALARS_H + i] > 0)
                    tcdc[i] = 100*1;
            }
        }
        for (int k = 0; k < NO_OF_SOLID_TRACERS; ++k)
            sprate[i] += fmax(ret_sink_velocity(k, 0, 0.001)*state_write_out -> tracer_densities[k*NO_OF_SCALARS + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i], 0);
        for (int k = NO_OF_SOLID_TRACERS; k < NO_OF_CONDENSATED_TRACERS; ++k)
            rprate[i] += fmax(ret_sink_velocity(k, 0, 0.001)*state_write_out -> tracer_densities[k*NO_OF_SCALARS + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i], 0);
    }
    
    // 10 m wind diagnostics
	double *wind_10_m_both_dir_array = malloc(2*min_no_of_output_steps*NO_OF_VECTORS_H*sizeof(double));
	double wind_10_m_downscale_factor = 0.667;
	double wind_tangential, wind_u_value, wind_v_value;
	int time_step_10_m_wind, h_index;
	double *wind_10_m_speed = malloc(min_no_of_output_steps*NO_OF_VECTORS_H*sizeof(double));
	double *wind_10_m_mean_u = malloc(NO_OF_VECTORS_H*sizeof(double));
	double *wind_10_m_mean_v = malloc(NO_OF_VECTORS_H*sizeof(double));
	for (int j = 0; j < min_no_of_output_steps*NO_OF_VECTORS_H; ++j)
	{
		time_step_10_m_wind = j/NO_OF_VECTORS_H;
		h_index = j - time_step_10_m_wind*NO_OF_VECTORS_H;
		wind_10_m_both_dir_array[2*j + 0] = wind_10_m_downscale_factor*wind_h_lowest_layer_array[j];
		wind_tangential = 0;
		for (int i = 0; i < 10; ++i)
			wind_tangential += grid -> trsk_modified_weights[10*h_index + i]*wind_h_lowest_layer_array[time_step_10_m_wind*NO_OF_VECTORS_H + grid -> trsk_modified_velocity_indices[10*h_index + i]];
		wind_10_m_both_dir_array[2*j + 1] = wind_10_m_downscale_factor*wind_tangential;
		wind_10_m_speed[j] = sqrt(pow(wind_10_m_both_dir_array[2*j + 0], 2) + pow(wind_10_m_both_dir_array[2*j + 1], 2));
		if (time_step_10_m_wind == 0)
		{
			wind_10_m_mean_u[h_index] = 1.0/min_no_of_output_steps*wind_10_m_both_dir_array[2*j + 0];
			wind_10_m_mean_v[h_index] = 1.0/min_no_of_output_steps*wind_10_m_both_dir_array[2*j + 1];
		}
		else
		{
			wind_10_m_mean_u[h_index] += 1.0/min_no_of_output_steps*wind_10_m_both_dir_array[2*j + 0];
			wind_10_m_mean_v[h_index] += 1.0/min_no_of_output_steps*wind_10_m_both_dir_array[2*j + 1];
		}
	}
	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
	{
	    passive_turn(wind_10_m_mean_u[i], wind_10_m_mean_v[i], -grid -> direction[i], &wind_u_value, &wind_v_value);
	    wind_10_m_mean_u[i] = wind_u_value;
	    wind_10_m_mean_v[i] = wind_v_value;
	}
	double standard_deviation;
	double gusts_parameter = 3;
	double *wind_10_m_gusts_speed = malloc(NO_OF_VECTORS_H*sizeof(double));
	double *vector_for_std_deviation = malloc(min_no_of_output_steps*sizeof(double));
	double wind_speed_10_m_mean;
	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		wind_speed_10_m_mean = 0;
		for (int j = 0; j < min_no_of_output_steps; ++j)
		{
			vector_for_std_deviation[j] = wind_10_m_speed[j*NO_OF_VECTORS_H + i];
			wind_speed_10_m_mean += 1.0/min_no_of_output_steps*wind_10_m_speed[j*NO_OF_VECTORS_H + i];
		}
		standard_deviation = calc_std_dev(vector_for_std_deviation, min_no_of_output_steps);
		if (t_write != t_init)
			wind_10_m_gusts_speed[i] = wind_speed_10_m_mean + gusts_parameter*standard_deviation;
		else
			wind_10_m_gusts_speed[i] = (1 + 0.2)*wind_speed_10_m_mean;
	}
	free(vector_for_std_deviation);
    
    // Diagnostics of quantities that are not surface-specific.    
	double *divv_h_all_layers = malloc(NO_OF_SCALARS*sizeof(double));
	if (write_out_divv_h == 1)
	{
		divv_h(state_write_out -> velocity_gas, divv_h_all_layers, grid);
	}
    Curl_field *rel_vort = calloc(1, sizeof(Curl_field));
	calc_rel_vort(state_write_out -> velocity_gas, *rel_vort, grid, dualgrid);
	
	// Diagnozing the u and v wind components at the vector points.
    double wind_0, wind_1;
	double *wind_u = malloc(NO_OF_H_VECTORS*sizeof(double));
	double *wind_v = malloc(NO_OF_H_VECTORS*sizeof(double));
	double *wind_w = malloc(NO_OF_V_VECTORS*sizeof(double));
	for (int i = 0; i < NO_OF_H_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
        wind_0 = state_write_out -> velocity_gas[h_index + layer_index*NO_OF_VECTORS_H + (layer_index + 1)*NO_OF_SCALARS_H];
        recov_hor_par_pri(state_write_out -> velocity_gas, layer_index, h_index, &wind_1, grid);
        passive_turn(wind_0, wind_1, -grid -> direction[h_index], &wind_u_value, &wind_v_value);
        wind_u[i] = wind_u_value;
        wind_v[i] = wind_v_value;
	}
	for (int i = 0; i < NO_OF_V_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_PER_LAYER;
		h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
		wind_w[i] = state_write_out -> velocity_gas[layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
	}
    Scalar_field *rh = calloc(1, sizeof(Scalar_field));
    Scalar_field *pressure = calloc(1, sizeof(Scalar_field));
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	(*rh)[i] = 100*rel_humidity(state_write_out -> temperature_gas[i], state_write_out -> tracer_densities[2*NO_OF_SCALARS + i]);
    	(*pressure)[i] = (state_write_out -> tracer_densities[2*NO_OF_SCALARS + i]*R_V + state_write_out -> density_dry[i]*R_D)*state_write_out -> temperature_gas[i];
    }
    
    int retval;
    if (io_config -> synop_output_switch == 1)
    {
    	double (*geopotential_height)[NO_OF_PRESSURE_LEVELS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_PRESSURE_LEVELS]));
    	double (*t_on_pressure_levels)[NO_OF_PRESSURE_LEVELS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_PRESSURE_LEVELS]));
    	double (*rh_on_pressure_levels)[NO_OF_PRESSURE_LEVELS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_PRESSURE_LEVELS]));
    	double (*u_on_pressure_levels)[NO_OF_PRESSURE_LEVELS] = malloc(sizeof(double[NO_OF_VECTORS_H][NO_OF_PRESSURE_LEVELS]));
    	double (*v_on_pressure_levels)[NO_OF_PRESSURE_LEVELS] = malloc(sizeof(double[NO_OF_VECTORS_H][NO_OF_PRESSURE_LEVELS]));
		// Netcdf output.
		if (io_config -> netcdf_output_switch == 1)
		{
			int OUTPUT_FILE_SYNOP_LENGTH = 300;
			char *OUTPUT_FILE_SYNOP_PRE = malloc((OUTPUT_FILE_SYNOP_LENGTH + 1)*sizeof(char));
			sprintf(OUTPUT_FILE_SYNOP_PRE, "output/%s/%s+%ds_synop.nc", RUN_ID, RUN_ID, (int) (t_write - t_init));
			OUTPUT_FILE_SYNOP_LENGTH = strlen(OUTPUT_FILE_SYNOP_PRE);
			free(OUTPUT_FILE_SYNOP_PRE);
			char *OUTPUT_FILE_SYNOP = malloc((OUTPUT_FILE_SYNOP_LENGTH + 1)*sizeof(char));
			sprintf(OUTPUT_FILE_SYNOP, "output/%s/%s+%ds_synop.nc", RUN_ID, RUN_ID, (int) (t_write - t_init));
			int ncid_synop, scalar_h_dimid, vector_h_dimid, level_dimid, geopot_height_id, temp_synop_id, rh_synop_id, wind_u_synop_id, wind_v_synop_id;
			
			if ((retval = nc_create(OUTPUT_FILE_SYNOP, NC_CLOBBER, &ncid_synop)))
				NCERR(retval);
			free(OUTPUT_FILE_SYNOP);
			if ((retval = nc_def_dim(ncid_synop, "scalar_index_h", NO_OF_SCALARS_H, &scalar_h_dimid)))
				NCERR(retval);
			if ((retval = nc_def_dim(ncid_synop, "vector_index_h", NO_OF_VECTORS_H, &vector_h_dimid)))
				NCERR(retval);
			if ((retval = nc_def_dim(ncid_synop, "level_index", NO_OF_PRESSURE_LEVELS, &level_dimid)))
				NCERR(retval);
			int dimids_synop_scalar[2];
			dimids_synop_scalar[0] = scalar_h_dimid;
			dimids_synop_scalar[1] = level_dimid;
			int dimids_synop_vector[2];
			dimids_synop_vector[0] = vector_h_dimid;
			dimids_synop_vector[1] = level_dimid;
			// Defining the variables.
			if ((retval = nc_def_var(ncid_synop, "geopotential_height", NC_DOUBLE, 2, dimids_synop_scalar, &geopot_height_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_synop, geopot_height_id, "units", strlen("gpm"), "gpm")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_synop, "temperature", NC_DOUBLE, 2, dimids_synop_scalar, &temp_synop_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_synop, temp_synop_id, "units", strlen("K"), "K")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_synop, "relative_humidity", NC_DOUBLE, 2, dimids_synop_scalar, &rh_synop_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_synop, rh_synop_id, "units", strlen("%"), "%")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_synop, "wind_u", NC_DOUBLE, 2, dimids_synop_vector, &wind_u_synop_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_synop, wind_u_synop_id, "units", strlen("m/s"), "m/s")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_synop, "wind_v", NC_DOUBLE, 2, dimids_synop_vector, &wind_v_synop_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_synop, wind_v_synop_id, "units", strlen("m/s"), "m/s")))
				NCERR(retval);
			if ((retval = nc_enddef(ncid_synop)))
				NCERR(retval);
			
			// Writing the arrays.
			if ((retval = nc_put_var_double(ncid_synop, geopot_height_id, &geopotential_height[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_synop, temp_synop_id, &t_on_pressure_levels[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_synop, rh_synop_id, &rh_on_pressure_levels[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_synop, wind_u_synop_id, &u_on_pressure_levels[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_synop, wind_v_synop_id, &v_on_pressure_levels[0][0])))
				NCERR(retval);
			
			// Closing the netcdf file.
			if ((retval = nc_close(ncid_synop)))
				NCERR(retval);
		}
    	free(geopotential_height);
    	free(t_on_pressure_levels);
    	free(rh_on_pressure_levels);
    	free(u_on_pressure_levels);
    	free(v_on_pressure_levels);
    }

    if (io_config -> aviation_output_switch == 1)
    {
    	double (*t_on_flight_levels)[NO_OF_FLIGHT_LEVELS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_FLIGHT_LEVELS]));
    	double (*rh_on_flight_levels)[NO_OF_FLIGHT_LEVELS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_FLIGHT_LEVELS]));
    	double (*u_on_flight_levels)[NO_OF_FLIGHT_LEVELS] = malloc(sizeof(double[NO_OF_VECTORS_H][NO_OF_FLIGHT_LEVELS]));
    	double (*v_on_flight_levels)[NO_OF_FLIGHT_LEVELS] = malloc(sizeof(double[NO_OF_VECTORS_H][NO_OF_FLIGHT_LEVELS]));
		// Netcdf output.
		if (io_config -> netcdf_output_switch == 1)
		{
			int OUTPUT_FILE_AVIATION_LENGTH = 300;
			char *OUTPUT_FILE_AVIATION_PRE = malloc((OUTPUT_FILE_AVIATION_LENGTH + 1)*sizeof(char));
			sprintf(OUTPUT_FILE_AVIATION_PRE, "output/%s/%s+%ds_aviation.nc", RUN_ID, RUN_ID, (int) (t_write - t_init));
			OUTPUT_FILE_AVIATION_LENGTH = strlen(OUTPUT_FILE_AVIATION_PRE);
			free(OUTPUT_FILE_AVIATION_PRE);
			char *OUTPUT_FILE_AVIATION = malloc((OUTPUT_FILE_AVIATION_LENGTH + 1)*sizeof(char));
			sprintf(OUTPUT_FILE_AVIATION, "output/%s/%s+%ds_aviation.nc", RUN_ID, RUN_ID, (int) (t_write - t_init));
			int ncid_aviation, scalar_h_dimid, vector_h_dimid, level_dimid, temp_aviation_id, rh_aviation_id, wind_u_aviation_id, wind_v_aviation_id;
			
			if ((retval = nc_create(OUTPUT_FILE_AVIATION, NC_CLOBBER, &ncid_aviation)))
				NCERR(retval);
			free(OUTPUT_FILE_AVIATION);
			if ((retval = nc_def_dim(ncid_aviation, "scalar_index_h", NO_OF_SCALARS_H, &scalar_h_dimid)))
				NCERR(retval);
			if ((retval = nc_def_dim(ncid_aviation, "vector_index_h", NO_OF_VECTORS_H, &vector_h_dimid)))
				NCERR(retval);
			if ((retval = nc_def_dim(ncid_aviation, "level_index", NO_OF_FLIGHT_LEVELS, &level_dimid)))
				NCERR(retval);
			int dimids_aviation_scalar[2];
			dimids_aviation_scalar[0] = scalar_h_dimid;
			dimids_aviation_scalar[1] = level_dimid;
			int dimids_aviation_vector[2];
			dimids_aviation_vector[0] = vector_h_dimid;
			dimids_aviation_vector[1] = level_dimid;
			
			// Defining the variables.
			if ((retval = nc_def_var(ncid_aviation, "temperature", NC_DOUBLE, 2, dimids_aviation_scalar, &temp_aviation_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_aviation, temp_aviation_id, "units", strlen("K"), "K")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_aviation, "relative_humidity", NC_DOUBLE, 2, dimids_aviation_scalar, &rh_aviation_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_aviation, rh_aviation_id, "units", strlen("%"), "%")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_aviation, "wind_u", NC_DOUBLE, 2, dimids_aviation_vector, &wind_u_aviation_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_aviation, wind_u_aviation_id, "units", strlen("m/s"), "m/s")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_aviation, "wind_v", NC_DOUBLE, 2, dimids_aviation_vector, &wind_v_aviation_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_aviation, wind_v_aviation_id, "units", strlen("m/s"), "m/s")))
				NCERR(retval);
			if ((retval = nc_enddef(ncid_aviation)))
				NCERR(retval);
			
			// Writing the arrays.
			if ((retval = nc_put_var_double(ncid_aviation, temp_aviation_id, &t_on_flight_levels[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_aviation, rh_aviation_id, &rh_on_flight_levels[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_aviation, wind_u_aviation_id, &u_on_flight_levels[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_aviation, wind_v_aviation_id, &v_on_flight_levels[0][0])))
				NCERR(retval);
			
			// Closing the netcdf file.
			if ((retval = nc_close(ncid_aviation)))
				NCERR(retval);
		}
    	free(t_on_flight_levels);
    	free(rh_on_flight_levels);
    	free(u_on_flight_levels);
    	free(v_on_flight_levels);
    }
    
    // Grib output.
    if (io_config -> grib_output_switch == 1)
    {
    	long unsigned length = 4;
		// Grib requires everything to be on horizontal levels.
		double *temperature_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *pressure_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *rh_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *wind_u_h = malloc(NO_OF_VECTORS_H*sizeof(double));
		double *wind_v_h = malloc(NO_OF_VECTORS_H*sizeof(double));
		double *rel_vort_h = malloc(NO_OF_VECTORS_H*sizeof(double));
		double *divv_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *wind_w_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		int OUTPUT_FILE_LENGTH = 300;
		char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
		sprintf(OUTPUT_FILE_PRE, "output/%s/%s+%ds.grb2", RUN_ID, RUN_ID, (int) (t_write - t_init));
		OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
		free(OUTPUT_FILE_PRE);
		char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
		sprintf(OUTPUT_FILE, "output/%s/%s+%ds.grb2", RUN_ID, RUN_ID, (int) (t_write - t_init));
		char *SAMPLE_FILENAME = "./input/grib_template.grb2";
		FILE *SAMPLE_FILE;
		int err = 0;
		if (t_init < 0)
			exit(1);
		FILE *OUT_GRIB;
		OUT_GRIB = fopen(OUTPUT_FILE, "w+");
		codes_handle *handle_temperature_h = NULL;
		codes_handle *handle_pressure_h = NULL;
		codes_handle *handle_wind_u_h = NULL;
		codes_handle *handle_wind_v_h = NULL;
		codes_handle *handle_wind_w_h = NULL;
		codes_handle *handle_wind_u_10m_mean = NULL;
		codes_handle *handle_wind_v_10m_mean = NULL;
		codes_handle *handle_mslp = NULL;
		codes_handle *handle_surface_p = NULL;
		codes_handle *handle_t2 = NULL;
		codes_handle *handle_tcdc = NULL;
		codes_handle *handle_rprate = NULL;
		codes_handle *handle_sprate = NULL;
		codes_handle *handle_wind_10m_gusts = NULL;
		codes_handle *handle_rel_vort = NULL;
		codes_handle *handle_rh = NULL;
		codes_handle *handle_cape = NULL;
		codes_handle *handle_divv_h = NULL;
		for (int i = 0; i < NO_OF_LAYERS; ++i)
		{
		    for (int j = 0; j < NO_OF_SCALARS_H; ++j)
		    {
		    	temperature_h[j] = state_write_out -> temperature_gas[i*NO_OF_SCALARS_H + j];
		    	pressure_h[j] = (*pressure)[i*NO_OF_SCALARS_H + j];
		    	rh_h[j] = (*rh)[i*NO_OF_SCALARS_H + j];
		    }
		    if (i == NO_OF_LAYERS - 1)
		    {
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_mslp = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
		        if ((retval = codes_set_long(handle_mslp, "discipline", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "centre", 255)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "significanceOfReferenceTime", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "productionStatusOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "typeOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "indicatorOfUnitOfTimeRange", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "stepUnits", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "dataDate", data_date)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "dataTime", data_time)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "forecastTime", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "stepRange", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "typeOfGeneratingProcess", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "parameterCategory", 3)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "parameterNumber", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "typeOfFirstFixedSurface", 102)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "scaledValueOfFirstFixedSurface", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "scaleFactorOfFirstFixedSurface", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_mslp, "level", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_double_array(handle_mslp, "values", mslp, NO_OF_SCALARS_H)))
		            ECCERR(retval);
		        if ((retval = codes_write_message(handle_mslp, OUTPUT_FILE, "a")))
		            ECCERR(retval);
				codes_handle_delete(handle_mslp);
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_surface_p = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
		        if ((retval = codes_write_message(handle_surface_p, OUTPUT_FILE, "a")))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "discipline", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "centre", 255)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "significanceOfReferenceTime", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "productionStatusOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "typeOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "indicatorOfUnitOfTimeRange", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "stepUnits", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "dataDate", data_date)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "dataTime", data_time)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "forecastTime", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "stepRange", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "typeOfGeneratingProcess", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "parameterCategory", 3)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "parameterNumber", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "typeOfFirstFixedSurface", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "scaledValueOfFirstFixedSurface", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "scaleFactorOfFirstFixedSurface", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_surface_p, "level", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_double_array(handle_surface_p, "values", surface_p, NO_OF_SCALARS_H)))
		            ECCERR(retval);
		        if ((retval = codes_write_message(handle_surface_p, OUTPUT_FILE, "a")))
		            ECCERR(retval);
				codes_handle_delete(handle_surface_p);
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_t2 = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
		        if ((retval = codes_set_long(handle_t2, "discipline", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "centre", 255)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "significanceOfReferenceTime", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "productionStatusOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "typeOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "indicatorOfUnitOfTimeRange", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "stepUnits", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "dataDate", data_date)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "dataTime", data_time)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "forecastTime", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "stepRange", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "typeOfGeneratingProcess", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "parameterCategory", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "parameterNumber", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "typeOfFirstFixedSurface", 103)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "scaledValueOfFirstFixedSurface", 2)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "scaleFactorOfFirstFixedSurface", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_t2, "level", 2)))
		            ECCERR(retval);
		        if ((retval = codes_set_double_array(handle_t2, "values", t2, NO_OF_SCALARS_H)))
		            ECCERR(retval);
		        if ((retval = codes_write_message(handle_t2, OUTPUT_FILE, "a")))
		            ECCERR(retval);
				codes_handle_delete(handle_t2);
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_rprate = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
		        if ((retval = codes_set_long(handle_rprate, "discipline", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "centre", 255)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "significanceOfReferenceTime", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "productionStatusOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "typeOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "indicatorOfUnitOfTimeRange", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "stepUnits", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "dataDate", data_date)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "dataTime", data_time)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "forecastTime", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "stepRange", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "typeOfGeneratingProcess", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "parameterCategory", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "parameterNumber", 65)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "typeOfFirstFixedSurface", 103)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "scaledValueOfFirstFixedSurface", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "scaleFactorOfFirstFixedSurface", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_rprate, "level", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_double_array(handle_rprate, "values", rprate, NO_OF_SCALARS_H)))
		            ECCERR(retval);
		        if ((retval = codes_write_message(handle_rprate, OUTPUT_FILE, "a")))
		            ECCERR(retval);
				codes_handle_delete(handle_rprate);
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_cape = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
		        if ((retval = codes_set_long(handle_cape, "discipline", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "centre", 255)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "significanceOfReferenceTime", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "productionStatusOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "typeOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "indicatorOfUnitOfTimeRange", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "stepUnits", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "dataDate", data_date)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "dataTime", data_time)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "forecastTime", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "stepRange", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "typeOfGeneratingProcess", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "parameterCategory", 7)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "parameterNumber", 6)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "typeOfFirstFixedSurface", 103)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "scaledValueOfFirstFixedSurface", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "scaleFactorOfFirstFixedSurface", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_cape, "level", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_double_array(handle_cape, "values", cape, NO_OF_SCALARS_H)))
		            ECCERR(retval);
		        if ((retval = codes_write_message(handle_cape, OUTPUT_FILE, "a")))
		            ECCERR(retval);
				codes_handle_delete(handle_cape);
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_sprate = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
		        if ((retval = codes_set_long(handle_sprate, "discipline", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "centre", 255)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "significanceOfReferenceTime", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "productionStatusOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "typeOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "indicatorOfUnitOfTimeRange", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "stepUnits", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "dataDate", data_date)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "dataTime", data_time)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "forecastTime", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "stepRange", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "typeOfGeneratingProcess", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "parameterCategory", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "parameterNumber", 66)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "typeOfFirstFixedSurface", 103)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "scaledValueOfFirstFixedSurface", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "scaleFactorOfFirstFixedSurface", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_sprate, "level", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_double_array(handle_sprate, "values", sprate, NO_OF_SCALARS_H)))
		            ECCERR(retval);
		        if ((retval = codes_write_message(handle_sprate, OUTPUT_FILE, "a")))
		            ECCERR(retval);
				codes_handle_delete(handle_sprate);
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_tcdc = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
		        if ((retval = codes_set_long(handle_tcdc, "discipline", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "centre", 255)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "significanceOfReferenceTime", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "productionStatusOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "typeOfProcessedData", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "indicatorOfUnitOfTimeRange", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "stepUnits", 13)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "dataDate", data_date)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "dataTime", data_time)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "forecastTime", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "stepRange", t_write - t_init)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "typeOfGeneratingProcess", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "parameterCategory", 6)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "parameterNumber", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "typeOfFirstFixedSurface", 103)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "scaledValueOfFirstFixedSurface", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "scaleFactorOfFirstFixedSurface", 1)))
		            ECCERR(retval);
		        if ((retval = codes_set_long(handle_tcdc, "level", 0)))
		            ECCERR(retval);
		        if ((retval = codes_set_string(handle_tcdc, "shortName", "tcc", &length)))
		            ECCERR(retval);
		        if ((retval = codes_set_double_array(handle_tcdc, "values", tcdc, NO_OF_SCALARS_H)))
		            ECCERR(retval);
		        if ((retval = codes_write_message(handle_tcdc, OUTPUT_FILE, "a")))
		            ECCERR(retval);
				codes_handle_delete(handle_tcdc);
		    }
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_temperature_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
		    if ((retval = codes_set_long(handle_temperature_h, "discipline", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "centre", 255)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "significanceOfReferenceTime", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "productionStatusOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "typeOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "indicatorOfUnitOfTimeRange", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "stepUnits", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "dataDate", data_date)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "dataTime", data_time)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "forecastTime", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "stepRange", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "typeOfGeneratingProcess", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "parameterCategory", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "parameterNumber", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "typeOfFirstFixedSurface", 26)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "scaledValueOfFirstFixedSurface", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_temperature_h, "level", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_double_array(handle_temperature_h, "values", temperature_h, NO_OF_SCALARS_H)))
		        ECCERR(retval);
		    if (i == 0)
		    {
		        if ((retval = codes_write_message(handle_temperature_h, OUTPUT_FILE, "w")))
		            ECCERR(retval);
		    }
		    else
		    {
		        if ((retval = codes_write_message(handle_temperature_h, OUTPUT_FILE, "a")))
		            ECCERR(retval);
		    }
			codes_handle_delete(handle_temperature_h);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_pressure_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
		    if ((retval = codes_set_long(handle_pressure_h, "discipline", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "centre", 255)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "significanceOfReferenceTime", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "productionStatusOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "typeOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "indicatorOfUnitOfTimeRange", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "stepUnits", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "dataDate", data_date)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "dataTime", data_time)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "forecastTime", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "stepRange", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "typeOfGeneratingProcess", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "parameterCategory", 3)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "parameterNumber", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "typeOfFirstFixedSurface", 26)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "scaledValueOfFirstFixedSurface", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_pressure_h, "level", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_double_array(handle_pressure_h, "values", pressure_h, NO_OF_SCALARS_H)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_pressure_h, OUTPUT_FILE, "a")))
		        ECCERR(retval);
			codes_handle_delete(handle_pressure_h);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_rh = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
		    if ((retval = codes_set_long(handle_rh, "discipline", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "centre", 255)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "significanceOfReferenceTime", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "productionStatusOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "typeOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "indicatorOfUnitOfTimeRange", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "stepUnits", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "dataDate", data_date)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "dataTime", data_time)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "forecastTime", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "stepRange", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "typeOfGeneratingProcess", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "parameterCategory", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "parameterNumber", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "typeOfFirstFixedSurface", 26)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "scaledValueOfFirstFixedSurface", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rh, "level", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_double_array(handle_rh, "values", rh_h, NO_OF_SCALARS_H)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_rh, OUTPUT_FILE, "a")))
		        ECCERR(retval);
			codes_handle_delete(handle_rh);
		    for (int j = 0; j < NO_OF_VECTORS_H; ++j)
		    {
		        wind_u_h[j] = wind_u[i*NO_OF_VECTORS_H + j];
		        wind_v_h[j] = wind_v[i*NO_OF_VECTORS_H + j];
		        rel_vort_h[j] = (*rel_vort)[NO_OF_VECTORS_H + i*2*NO_OF_VECTORS_H + j];
		    }
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_wind_u_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
		    if ((retval = codes_set_long(handle_wind_u_h, "discipline", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "centre", 255)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "significanceOfReferenceTime", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "productionStatusOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "typeOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "indicatorOfUnitOfTimeRange", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "stepUnits", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "dataDate", data_date)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "dataTime", data_time)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "forecastTime", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "stepRange", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "typeOfGeneratingProcess", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "parameterCategory", 2)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "parameterNumber", 2)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "typeOfFirstFixedSurface", 26)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "scaledValueOfFirstFixedSurface", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_u_h, "level", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_double_array(handle_wind_u_h, "values", wind_u_h, NO_OF_VECTORS_H)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_wind_u_h, OUTPUT_FILE, "a")))
		        ECCERR(retval);
			codes_handle_delete(handle_wind_u_h);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_wind_v_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
		    if ((retval = codes_set_long(handle_wind_v_h, "discipline", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "centre", 255)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "significanceOfReferenceTime", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "productionStatusOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "typeOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "indicatorOfUnitOfTimeRange", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "stepUnits", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "dataDate", data_date)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "dataTime", data_time)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "forecastTime", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "stepRange", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "typeOfGeneratingProcess", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "parameterCategory", 2)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "parameterNumber", 3)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "typeOfFirstFixedSurface", 26)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "scaledValueOfFirstFixedSurface", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_v_h, "level", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_double_array(handle_wind_v_h, "values", wind_v_h, NO_OF_VECTORS_H)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_wind_v_h, OUTPUT_FILE, "a")))
		        ECCERR(retval);
			codes_handle_delete(handle_wind_v_h);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_rel_vort = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
		    if ((retval = codes_set_long(handle_rel_vort, "discipline", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "centre", 255)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "significanceOfReferenceTime", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "productionStatusOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "typeOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "indicatorOfUnitOfTimeRange", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "stepUnits", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "dataDate", data_date)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "dataTime", data_time)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "forecastTime", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "stepRange", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "typeOfGeneratingProcess", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "parameterCategory", 2)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "parameterNumber", 12)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "typeOfFirstFixedSurface", 26)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "scaledValueOfFirstFixedSurface", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rel_vort, "level", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_double_array(handle_rel_vort, "values", rel_vort_h, NO_OF_VECTORS_H)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_rel_vort, OUTPUT_FILE, "a")))
		        ECCERR(retval);
			codes_handle_delete(handle_rel_vort);
			if (write_out_divv_h == 1)
			{	
				for (int j = 0; j < NO_OF_SCALARS_H; ++j)
				{
					divv_h[j] = divv_h_all_layers[i*NO_OF_SCALARS_H + j];
				}
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_divv_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
				if ((retval = codes_set_long(handle_divv_h, "discipline", 0)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "centre", 255)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "significanceOfReferenceTime", 1)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "productionStatusOfProcessedData", 1)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "typeOfProcessedData", 1)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "indicatorOfUnitOfTimeRange", 13)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "stepUnits", 13)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "dataDate", data_date)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "dataTime", data_time)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "forecastTime", t_write - t_init)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "stepRange", t_write - t_init)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "typeOfGeneratingProcess", 1)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "parameterCategory", 2)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "parameterNumber", 13)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "typeOfFirstFixedSurface", 26)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "scaledValueOfFirstFixedSurface", i)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "scaleFactorOfFirstFixedSurface", 1)))
				    ECCERR(retval);
				if ((retval = codes_set_long(handle_divv_h, "level", i)))
				    ECCERR(retval);
				if ((retval = codes_set_double_array(handle_divv_h, "values", divv_h, NO_OF_SCALARS_H)))
				    ECCERR(retval);
				if ((retval = codes_write_message(handle_divv_h, OUTPUT_FILE, "a")))
				    ECCERR(retval);
				codes_handle_delete(handle_divv_h);
			}
		}
		SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_wind_u_10m_mean = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		fclose(SAMPLE_FILE);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "discipline", 0)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "centre", 255)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "significanceOfReferenceTime", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "productionStatusOfProcessedData", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "typeOfProcessedData", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "indicatorOfUnitOfTimeRange", 13)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "stepUnits", 13)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "dataDate", data_date)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "dataTime", data_time)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "forecastTime", t_write - t_init)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "stepRange", t_write - t_init)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "typeOfGeneratingProcess", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "parameterCategory", 2)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "parameterNumber", 2)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "typeOfFirstFixedSurface", 103)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "scaledValueOfFirstFixedSurface", 10)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "scaleFactorOfFirstFixedSurface", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_u_10m_mean, "level", 10)))
			ECCERR(retval);
		if ((retval = codes_set_double_array(handle_wind_u_10m_mean, "values", wind_10_m_mean_u, NO_OF_VECTORS_H)))
			ECCERR(retval);
		if ((retval = codes_write_message(handle_wind_u_10m_mean, OUTPUT_FILE, "a")))
			ECCERR(retval);
		codes_handle_delete(handle_wind_u_10m_mean);
		free(wind_10_m_mean_u);
		SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_wind_v_10m_mean = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		fclose(SAMPLE_FILE);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "discipline", 0)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "centre", 255)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "significanceOfReferenceTime", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "productionStatusOfProcessedData", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "typeOfProcessedData", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "indicatorOfUnitOfTimeRange", 13)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "stepUnits", 13)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "dataDate", data_date)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "dataTime", data_time)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "forecastTime", t_write - t_init)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "stepRange", t_write - t_init)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "typeOfGeneratingProcess", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "parameterCategory", 2)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "parameterNumber", 3)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "typeOfFirstFixedSurface", 103)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "scaledValueOfFirstFixedSurface", 10)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "scaleFactorOfFirstFixedSurface", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_v_10m_mean, "level", 10)))
			ECCERR(retval);
		if ((retval = codes_set_double_array(handle_wind_v_10m_mean, "values", wind_10_m_mean_v, NO_OF_VECTORS_H)))
			ECCERR(retval);
		if ((retval = codes_write_message(handle_wind_v_10m_mean, OUTPUT_FILE, "a")))
			ECCERR(retval);
		free(wind_10_m_mean_v);
		codes_handle_delete(handle_wind_v_10m_mean);
		SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_wind_10m_gusts = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		fclose(SAMPLE_FILE);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "discipline", 0)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "centre", 255)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "significanceOfReferenceTime", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "productionStatusOfProcessedData", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "typeOfProcessedData", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "indicatorOfUnitOfTimeRange", 13)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "stepUnits", 13)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "dataDate", data_date)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "dataTime", data_time)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "forecastTime", t_write - t_init)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "stepRange", t_write - t_init)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "typeOfGeneratingProcess", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "parameterCategory", 2)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "parameterNumber", 22)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "typeOfFirstFixedSurface", 103)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "scaledValueOfFirstFixedSurface", 10)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "scaleFactorOfFirstFixedSurface", 1)))
			ECCERR(retval);
		if ((retval = codes_set_long(handle_wind_10m_gusts, "level", 10)))
			ECCERR(retval);
		if ((retval = codes_set_double_array(handle_wind_10m_gusts, "values", wind_10_m_gusts_speed, NO_OF_VECTORS_H)))
			ECCERR(retval);
		if ((retval = codes_write_message(handle_wind_10m_gusts, OUTPUT_FILE, "a")))
			ECCERR(retval);
		codes_handle_delete(handle_wind_10m_gusts);
		free(wind_u_h);
		free(wind_v_h);
		free(rel_vort_h);
		free(divv_h);
		free(temperature_h);
		free(pressure_h);
		free(rh_h);
		SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_wind_w_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		fclose(SAMPLE_FILE);
		for (int i = 0; i < NO_OF_LEVELS; ++i)
		{
		    for (int j = 0; j < NO_OF_SCALARS_H; j++)
		    {
		        wind_w_h[j] = state_write_out -> velocity_gas[j + i*NO_OF_VECTORS_PER_LAYER];
		    }
		    if ((retval = codes_set_long(handle_wind_w_h, "discipline", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "centre", 255)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "significanceOfReferenceTime", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "productionStatusOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "typeOfProcessedData", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "indicatorOfUnitOfTimeRange", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "stepUnits", 13)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "dataDate", data_date)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "dataTime", data_time)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "forecastTime", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "stepRange", t_write - t_init)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "typeOfGeneratingProcess", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "parameterCategory", 2)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "parameterNumber", 9)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "typeOfFirstFixedSurface", 26)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "scaledValueOfFirstFixedSurface", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_wind_w_h, "level", i)))
		        ECCERR(retval);
		    if ((retval = codes_set_double_array(handle_wind_w_h, "values", wind_w_h, NO_OF_SCALARS_H)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_wind_w_h, OUTPUT_FILE, "a")))
		        ECCERR(retval);
		}
		free(OUTPUT_FILE);
		codes_handle_delete(handle_wind_w_h);
		free(wind_w_h);
		fclose(OUT_GRIB);
    }
    
    // Netcdf output.
    if (io_config -> netcdf_output_switch == 1)
    {
		int OUTPUT_FILE_LENGTH = 300;
		char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
		sprintf(OUTPUT_FILE_PRE, "output/%s/%s+%ds.nc", RUN_ID, RUN_ID, (int) (t_write - t_init));
		OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
		free(OUTPUT_FILE_PRE);
		char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
		sprintf(OUTPUT_FILE, "output/%s/%s+%ds.nc", RUN_ID, RUN_ID, (int) (t_write - t_init));
		int scalar_dimid, scalar_h_dimid, vector_h_dimid, vector_v_dimid, temp_id, density_dry_id, pressure_id, wind_u_id, wind_v_id, wind_w_id, rh_id, mslp_id, ncid, retval, surface_p_id, rprate_id, sprate_id, cape_id, tcdc_id, t2_id, divv_h_all_layers_id, rel_vort_id, curl_field_dimid;
		
		if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
		    NCERR(retval);
		free(OUTPUT_FILE);
		if ((retval = nc_def_dim(ncid, "scalar_index", NO_OF_SCALARS, &scalar_dimid)))
		    NCERR(retval);
		if ((retval = nc_def_dim(ncid, "scalar_index_h", NO_OF_SCALARS_H, &scalar_h_dimid)))
		    NCERR(retval);
		if ((retval = nc_def_dim(ncid, "vector_index_h", NO_OF_H_VECTORS, &vector_h_dimid)))
		    NCERR(retval);
		if ((retval = nc_def_dim(ncid, "vector_index_v", NO_OF_V_VECTORS, &vector_v_dimid)))
		    NCERR(retval);
		if ((retval = nc_def_dim(ncid, "curl_point_index", NO_OF_LAYERS*2*NO_OF_VECTORS_H + NO_OF_VECTORS_H, &curl_field_dimid)))
		    NCERR(retval);
		
		// Defining the variables.
		if ((retval = nc_def_var(ncid, "temperature_gas", NC_DOUBLE, 1, &scalar_dimid, &temp_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, temp_id, "units", strlen("K"), "K")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "density_dry", NC_DOUBLE, 1, &scalar_dimid, &density_dry_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, density_dry_id, "units", strlen("kg/m^3"), "kg/m^3")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "pressure", NC_DOUBLE, 1, &scalar_dimid, &pressure_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, pressure_id, "units", strlen("Pa"), "Pa")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "wind_u", NC_DOUBLE, 1, &vector_h_dimid, &wind_u_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, wind_u_id, "units", strlen("m/s"), "m/s")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "wind_v", NC_DOUBLE, 1, &vector_h_dimid, &wind_v_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, wind_v_id, "units", strlen("m/s"), "m/s")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "wind_w", NC_DOUBLE, 1, &vector_v_dimid, &wind_w_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, wind_w_id, "units", strlen("m/s"), "m/s")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "rh", NC_DOUBLE, 1, &scalar_dimid, &rh_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, rh_id, "units", strlen("%"), "%")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "rel_vort", NC_DOUBLE, 1, &curl_field_dimid, &rel_vort_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, rel_vort_id, "units", strlen("1/s"), "1/s")))
		    NCERR(retval);
		if (write_out_divv_h == 1)
		{
			if ((retval = nc_def_var(ncid, "divv_h_all_layers", NC_DOUBLE, 1, &scalar_dimid, &divv_h_all_layers_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid, divv_h_all_layers_id, "units", strlen("1/s"), "1/s")))
				NCERR(retval);
		}
		if ((retval = nc_def_var(ncid, "mslp", NC_DOUBLE, 1, &scalar_h_dimid, &mslp_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, mslp_id, "units", strlen("Pa"), "Pa")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "surface_p", NC_DOUBLE, 1, &scalar_h_dimid, &surface_p_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, surface_p_id, "units", strlen("Pa"), "Pa")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "t2", NC_DOUBLE, 1, &scalar_h_dimid, &t2_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, t2_id, "units", strlen("K"), "K")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "tcdc", NC_DOUBLE, 1, &scalar_h_dimid, &tcdc_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, tcdc_id, "units", strlen("%"), "%")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "rprate", NC_DOUBLE, 1, &scalar_h_dimid, &rprate_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, rprate_id, "units", strlen("kg/(m^2s)"), "kg/(m^2s)")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "sprate", NC_DOUBLE, 1, &scalar_h_dimid, &sprate_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, sprate_id, "units", strlen("kg/(m^2s)"), "kg/(m^2s)")))
		    NCERR(retval);
		if ((retval = nc_def_var(ncid, "cape", NC_DOUBLE, 1, &scalar_h_dimid, &cape_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, cape_id, "units", strlen("J/kg"), "J/kg")))
		    NCERR(retval);
		if ((retval = nc_enddef(ncid)))
		    NCERR(retval);
		
		// Writing the 3D arrays.
		if ((retval = nc_put_var_double(ncid, temp_id, &state_write_out -> temperature_gas[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, density_dry_id, &state_write_out -> density_dry[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, pressure_id, &(*rel_vort)[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, wind_u_id, &wind_u[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, wind_v_id, &wind_v[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, wind_w_id, &wind_w[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, rh_id, &(*rh)[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, rel_vort_id, &(*rel_vort)[0])))
		    NCERR(retval);
		if (write_out_divv_h == 1)
		{
			if ((retval = nc_put_var_double(ncid, divv_h_all_layers_id, &divv_h_all_layers[0])))
				NCERR(retval);
	    }
		    
		// Writing out the surface specific variables.
		if ((retval = nc_put_var_double(ncid, mslp_id, &mslp[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, surface_p_id, &surface_p[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, t2_id, &t2[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, tcdc_id, &tcdc[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, rprate_id, &rprate[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, sprate_id, &sprate[0])))
		    NCERR(retval);
		if ((retval = nc_put_var_double(ncid, cape_id, &cape[0])))
		    NCERR(retval);
		
		// Closing the netcdf file.
		if ((retval = nc_close(ncid)))
			NCERR(retval);
	}
	free(divv_h_all_layers);
	free(wind_10_m_speed);
	free(wind_10_m_gusts_speed);
	free(wind_10_m_both_dir_array);
	free(t2);
	free(mslp);
	free(surface_p);
	free(rprate);
	free(sprate);
	free(tcdc);
	free(pot_temperature);
	free(rel_vort);
	free(cape);
	free(wind_u);
	free(wind_v);
	free(wind_w);
	free(rh);
	free(pressure);
	return 0;
}

int write_out_integral(State *state_write_out, double t_write, char RUN_ID[], Grid *grid, Dualgrid *dualgrid, int integral_id)
{
	/*
	integral_id:
	0: dry mass
	1: entropy
	2: energy
	*/
    double global_integral = 0;
    FILE *global_integral_file;
    int INTEGRAL_FILE_LENGTH = 200;
    char *INTEGRAL_FILE_PRE = malloc((INTEGRAL_FILE_LENGTH + 1)*sizeof(char));
    if (integral_id == 0)
   		sprintf(INTEGRAL_FILE_PRE, "output/%s/%s", RUN_ID, "dry_mass");
    if (integral_id == 1)
   		sprintf(INTEGRAL_FILE_PRE, "output/%s/%s", RUN_ID, "entropy");
    if (integral_id == 2)
   		sprintf(INTEGRAL_FILE_PRE, "output/%s/%s", RUN_ID, "energy");
    if (integral_id == 2)
   		sprintf(INTEGRAL_FILE_PRE, "output/%s/%s", RUN_ID, "energy");
    if (integral_id == 3)
   		sprintf(INTEGRAL_FILE_PRE, "output/%s/%s", RUN_ID, "linearized_entropy");
    INTEGRAL_FILE_LENGTH = strlen(INTEGRAL_FILE_PRE);
    char *INTEGRAL_FILE = malloc((INTEGRAL_FILE_LENGTH + 1)*sizeof(char));
    sprintf(INTEGRAL_FILE, "%s", INTEGRAL_FILE_PRE);
    free(INTEGRAL_FILE_PRE);
    if (integral_id == 0)
    {
    	global_integral_file = fopen(INTEGRAL_FILE, "a");
    	global_scalar_integrator(state_write_out -> density_dry, grid, &global_integral);
    	fprintf(global_integral_file, "%lf\t%lf\n", t_write, global_integral);
    	fclose(global_integral_file);
    }
    if (integral_id == 1)
    {
    	global_integral_file = fopen(INTEGRAL_FILE, "a");
    	global_scalar_integrator(state_write_out -> entropy_density_dry, grid, &global_integral);
    	fprintf(global_integral_file, "%lf\t%lf\n", t_write, global_integral);
    	fclose(global_integral_file);
    }
    if (integral_id == 2)
    {
    	double kinetic_integral, potential_integral, internal_integral;
    	global_integral_file = fopen(INTEGRAL_FILE, "a");
    	Scalar_field *e_kin_density = malloc(sizeof(Scalar_field));
    	kinetic_energy(state_write_out -> velocity_gas, *e_kin_density, grid, 1);
    	scalar_times_scalar(state_write_out -> density_dry, *e_kin_density, *e_kin_density);
    	global_scalar_integrator(*e_kin_density, grid, &kinetic_integral);
    	free(e_kin_density);
    	Scalar_field *pot_energy_density = malloc(sizeof(Scalar_field));
    	scalar_times_scalar(state_write_out -> density_dry, grid -> gravity_potential, *pot_energy_density);
    	global_scalar_integrator(*pot_energy_density, grid, &potential_integral);
    	free(pot_energy_density);
    	Scalar_field *int_energy_density = malloc(sizeof(Scalar_field));
    	scalar_times_scalar(state_write_out -> density_dry, state_write_out -> temperature_gas, *int_energy_density);
    	global_scalar_integrator(*int_energy_density, grid, &internal_integral);
    	fprintf(global_integral_file, "%lf\t%lf\t%lf\t%lf\n", t_write, kinetic_integral, potential_integral, C_D_V*internal_integral);
    	free(int_energy_density);
    	fclose(global_integral_file);
    }
    if (integral_id == 3)
    {
    	global_integral_file = fopen(INTEGRAL_FILE, "a");
    	Scalar_field *pot_temp = malloc(sizeof(Scalar_field));
		pot_temp_diagnostics(state_write_out, *pot_temp);
    	Scalar_field *linear_entropy = malloc(sizeof(Scalar_field));
    	for (int i = 0; i < NO_OF_SCALARS; ++i)
    	{
    		(*linear_entropy)[i] = C_D_P*state_write_out -> density_dry[i]*(*pot_temp)[i];
    	}
    	global_scalar_integrator(*linear_entropy, grid, &global_integral);
    	fprintf(global_integral_file, "%lf\t%lf\n", t_write, global_integral);
    	free(linear_entropy);
    	free(pot_temp);
    	fclose(global_integral_file);
    }
    free(INTEGRAL_FILE);
	return 0;
}

double calc_std_dev(double vector_for_std_deviation[], int no_of_values)
{
	double mean = 0;
	for (int i = 0; i < no_of_values; ++i)
		mean += 1.0/no_of_values*vector_for_std_deviation[i];
	double result = 0;
	for (int i = 0; i < no_of_values; ++i)
		result += pow(vector_for_std_deviation[i] - mean, 2);
	result = 1/sqrt(no_of_values)*sqrt(result);
	return result;
}

