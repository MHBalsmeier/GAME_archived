/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
In this file, the initial state of the simulation is read in from a netcdf file.
*/

#include <stdio.h>
#include <stdlib.h>
#include "../enum_and_typedefs.h"
#include "io.h"
#include "../thermodynamics/thermodynamics.h"
#include "../settings.h"
#include <netcdf.h>
#include "atmostracers.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

int set_init_data(char FILE_NAME[], State *init_state, Grid* grid)
{
    double *density_dry = malloc(NO_OF_SCALARS*sizeof(double));
    double *temperature_gas = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *water_vapour_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_temperature = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_temperature = malloc(NO_OF_SCALARS*sizeof(double));
    int retval, ncid;
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        NCERR(retval);
    int temperature_gas_id, density_dry_id, wind_id, density_vapour_id, density_liquid_id, density_solid_id, temperature_liquid_id, temperature_solid_id, stretching_parameter_id;
    double stretching_parameter;
    if ((retval = nc_inq_varid(ncid, "density_dry", &density_dry_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature_gas",&temperature_gas_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "wind", &wind_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_solid", &density_solid_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_liquid", &density_liquid_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_vapour", &density_vapour_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature_liquid", &temperature_liquid_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature_solid", &temperature_solid_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "stretching_parameter", &stretching_parameter_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, stretching_parameter_id, &stretching_parameter)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, temperature_gas_id, &temperature_gas[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, density_dry_id, &density_dry[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, wind_id, &wind[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, density_vapour_id, &water_vapour_density[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, density_liquid_id, &liquid_water_density [0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, density_solid_id, &solid_water_density[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, temperature_liquid_id, &liquid_water_temperature[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, temperature_solid_id, &solid_water_temperature[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
        NCERR(retval);
    
    // checking wether the stretching parameters of the grid used for creating the input file and the grid file read in do conform
    if (grid -> stretching_parameter != stretching_parameter)
    {
    	printf("Stretching parameters of grid and input file do not conform.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
    
    // resricting the maximum relative humidity to 100 %
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
		if (rel_humidity(water_vapour_density[i], temperature_gas[i]) > 1)
		{
			water_vapour_density[i] = water_vapour_density[i]/rel_humidity(water_vapour_density[i], temperature_gas[i]);
		}
    }
    
    int layer_index, h_index;
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		if (NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_scalar[h_index])
		{
		    if (NO_OF_CONSTITUENTS >= 4)
		    {
				init_state -> rho[i] = solid_water_density[i];
				if (init_state -> rho[i] < 0)
				{
					printf("Error: init_state -> rho < 0 somewhere.\n");
					printf("Aborting.\n");
					exit(1);
				}
				init_state -> rho[NO_OF_SCALARS + i] = liquid_water_density[i];
				if (init_state -> rho[NO_OF_SCALARS + i] < 0)
				{
					printf("Error: init_state -> rho < 0 somewhere.\n");
					printf("Aborting.\n");
					exit(1);
				}
		    }
		    init_state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i] = density_dry[i];
		    if (init_state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i] < 0)
		    {
		    	printf("Error: init_state -> rho < 0 somewhere.\n");
		    	printf("Aborting.\n");
		    	exit(1);
		    }
		    if (NO_OF_CONSTITUENTS >= 4)
		    {
				init_state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i] = water_vapour_density[i];
				if (init_state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i] < 0)
				{
					printf("Error: init_state -> rho < 0 somewhere.\n");
					printf("Aborting.\n");
					exit(1);
				}
		    }
	    }
	}
	
    double pressure, pot_temp, specific_entropy;
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		if (NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_scalar[h_index])
		{
			for (int j = 0; j < NO_OF_CONSTITUENTS; ++j)
			{
				if (j == NO_OF_CONDENSED_CONSTITUENTS - 2)
				{
					init_state -> condensed_density_temperatures[i] = solid_water_density[i]*solid_water_temperature[i];
				}
				if (j == NO_OF_CONDENSED_CONSTITUENTS - 1)
				{
					init_state -> condensed_density_temperatures[NO_OF_SCALARS + i] = liquid_water_density[i]*liquid_water_temperature[i];
				}
				if (j >= NO_OF_CONDENSED_CONSTITUENTS)
				{
					if (init_state -> rho[j*NO_OF_SCALARS + i] == 0)
					{
						init_state -> rhotheta[(j - NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + i] = 0;
					}
					else
					{
						pressure = init_state -> rho[j*NO_OF_SCALARS + i]*specific_gas_constants(j - NO_OF_CONDENSED_CONSTITUENTS)*temperature_gas[i];
						pot_temp = temperature_gas[i]*pow(P_0/pressure, specific_gas_constants(j - NO_OF_CONDENSED_CONSTITUENTS)/spec_heat_capacities_p_gas(j - NO_OF_CONDENSED_CONSTITUENTS));
						specific_entropy = spec_heat_capacities_p_gas(j - NO_OF_CONDENSED_CONSTITUENTS)*log(pot_temp);
						init_state -> rhotheta[(j - NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + i]
						= init_state -> rho[j*NO_OF_SCALARS + i]*specific_entropy;
					}
				}
			}
	    }
	}
	
	for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	layer_index = i/NO_OF_VECTORS_PER_LAYER;
    	h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
    	// the horizontal case
		if (h_index >= NO_OF_SCALARS_H
		// check for shading
		&& NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_vector[h_index - NO_OF_SCALARS_H])
		{
        	init_state -> velocity_gas[i] = wind[i];
    	}
    	// the vertical case
    	if (h_index < NO_OF_SCALARS_H
    	// check for shading
    	&& NO_OF_LAYERS - layer_index > grid -> no_of_shaded_points_scalar[h_index])
        {
        	init_state -> velocity_gas[i] = wind[i];
    	}
    }
    free(density_dry);
    free(wind);
    free(water_vapour_density);
    free(liquid_water_density);
    free(solid_water_density);
    free(liquid_water_temperature);
    free(solid_water_temperature);
    return 0;
}






