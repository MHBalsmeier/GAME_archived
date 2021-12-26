/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, the initial state of the simulation is read in from a netcdf file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include <atmostracers.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../thermodynamics.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

int set_init_data(char FILE_NAME[], State *init_state, Grid* grid, Soil *soil)
{
	/*
	This function sets the initial state of the model atmosphere.
	*/
    double *temperatures = malloc((NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS*sizeof(double));
    double *sst = malloc(NO_OF_SCALARS_H*sizeof(double));
    int retval, ncid;
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        NCERR(retval);
    int densities_id, temperatures_id, wind_id, sst_id, stretching_parameter_id, sst_avail;
    double stretching_parameter;
    if ((retval = nc_inq_varid(ncid, "densities", &densities_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperatures", &temperatures_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "wind", &wind_id)))
        NCERR(retval);
    // figuring out if the netcdf file contains SST
    sst_avail = 0;
    if (nc_inq_varid(ncid, "sst", &sst_id) == 0)
    {
    	sst_avail = 1;
    }
    if ((retval = nc_inq_varid(ncid, "stretching_parameter", &stretching_parameter_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, densities_id, &init_state -> rho[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, temperatures_id, &temperatures[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, wind_id, &init_state -> wind[0])))
        NCERR(retval);
    // reading the SST data if it is present in the netcdf file
    if (sst_avail == 1)
    {
		if ((retval = nc_get_var_double(ncid, sst_id, &sst[0])))
		    NCERR(retval);
    }
    if ((retval = nc_get_var_double(ncid, stretching_parameter_id, &stretching_parameter)))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
        NCERR(retval);
    
	// resricting the maximum relative humidity to 100 %
    if (NO_OF_CONDENSED_CONSTITUENTS == 4)
    {
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			if (rel_humidity(init_state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i], temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]) > 1)
			{
				init_state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i] = init_state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i]
				/rel_humidity(init_state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i], temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]);
			}
		}
    }
	
    // determining the temperature densities of the condensates
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS; ++i)
	{
		init_state -> condensed_density_temperatures[i] = init_state -> rho[i]*temperatures[i];
	}
	
	// diagnostic thermodynamical quantities
	double pressure, pot_temp;
	#pragma omp parallel for private(pressure, pot_temp)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		pressure = init_state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*specific_gas_constants(0)*temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i];
		pot_temp = temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*pow(P_0/pressure, specific_gas_constants(0)/spec_heat_capacities_p_gas(0));
		init_state -> rhotheta[i] = init_state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*pot_temp;
		// calculating the potential temperature perturbation
		init_state -> theta_pert[i] = pot_temp - grid -> theta_bg[i];
		// calculating the Exner pressure perturbation
		init_state -> exner_pert[i] = temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]/(grid -> theta_bg[i] + init_state -> theta_pert[i]) - grid -> exner_bg[i];
	}
	
	int soil_index;
	double z_soil, t_sfc;
	// setting the soil temperature
	#pragma omp parallel for private(soil_index, z_soil, t_sfc)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// temperature at the surface
		if (grid -> is_land[i] == 1 || (grid -> is_land[i] == 0 && sst_avail == 0))
		{
			t_sfc = temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + NO_OF_SCALARS - NO_OF_SCALARS_H + i];
		}
		else
		{
			t_sfc = sst[i];
		}
		
		// loop over all soil layers
		for (int soil_layer_index = 0; soil_layer_index < NO_OF_SOIL_LAYERS; ++soil_layer_index)
		{
			// index of this soil grid point
			soil_index = i + soil_layer_index*NO_OF_SCALARS_H;
			z_soil = grid -> z_t_const/NO_OF_SOIL_LAYERS*(0.5 + soil_layer_index);
			soil -> temperature[soil_index] = t_sfc + (grid -> t_const_soil - t_sfc)*z_soil/grid -> z_t_const;
		}
	}
	
    // checks
    // checking for negative densities
    # pragma omp parallel for
	for (int i = 0; i < NO_OF_CONSTITUENTS*NO_OF_SCALARS; ++i)
	{
		if (init_state -> rho[i] < 0)
	    {
			printf("Negative density found.\n");
			printf("Aborting.\n");
			exit(1);
    	}
	}
	// checking wether the stretching parameters of the grid used for creating the input file and the grid file read in conform
    if (grid -> stretching_parameter != stretching_parameter)
    {
    	printf("Stretching parameters of grid and input file do not conform.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
    
    free(sst);
    free(temperatures);
    return 0;
}






