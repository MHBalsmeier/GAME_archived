/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include <stdio.h>
#include <stdlib.h>
#include "../enum_and_typedefs.h"
#include "io.h"
#include "../diagnostics/diagnostics.h"
#include "../settings.h"
#include <netcdf.h>
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

int set_init_data(char FILE_NAME[], State *init_state)
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
    int temperature_gas_id, density_dry_id, wind_id, density_vapour_id, density_liquid_id, density_solid_id, temperature_liquid_id, temperature_solid_id;
    if ((retval = nc_inq_varid(ncid, "density_dry", &density_dry_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature_gas",&temperature_gas_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "wind", &wind_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_vapour", &density_vapour_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_liquid", &density_liquid_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_solid", &density_solid_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature_liquid", &temperature_liquid_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature_solid", &temperature_solid_id)))
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
    
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        init_state -> temperature_gas[i] = temperature_gas[i];
    }
    
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
        init_state -> mass_densities[i] = solid_water_density[i];
        if (init_state -> mass_densities[i] < 0)
        {
        	printf("Error: init_state -> mass_densities < 0 somewhere.\n");
        	printf("Aborting.\n");
        	exit(1);
        }
        init_state -> mass_densities[NO_OF_SCALARS + i] = liquid_water_density[i];
        if (init_state -> mass_densities[NO_OF_SCALARS + i] < 0)
        {
        	printf("Error: init_state -> mass_densities < 0 somewhere.\n");
        	printf("Aborting.\n");
        	exit(1);
        }
        init_state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i] = density_dry[i];
        if (init_state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i] < 0)
        {
        	printf("Error: init_state -> mass_densities < 0 somewhere.\n");
        	printf("Aborting.\n");
        	exit(1);
        }
        init_state -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i] = water_vapour_density[i];
        if (init_state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i] < 0)
        {
        	printf("Error: init_state -> mass_densities < 0 somewhere.\n");
        	printf("Aborting.\n");
        	exit(1);
        }
	}
	
    double pressure, pot_temp, specific_entropy;
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		for (int j = 0; j < NO_OF_CONSTITUENTS; ++j)
		{
			if (j < NO_OF_CONDENSED_CONSTITUENTS)
			{
				init_state -> entropy_densities[i] = 0;
				init_state -> entropy_densities[NO_OF_SCALARS + i] = 0;
			}
			else
			{
				if (init_state -> mass_densities[j*NO_OF_SCALARS + i] == 0)
				{
					init_state -> entropy_densities[j*NO_OF_SCALARS + i] = 0;
				}
				else
				{
					pressure = init_state -> mass_densities[j*NO_OF_SCALARS + i]*specific_gas_constants(j - NO_OF_CONDENSED_CONSTITUENTS)*temperature_gas[i];
					pot_temp = temperature_gas[i]*pow(P_0/pressure, specific_gas_constants(j - NO_OF_CONDENSED_CONSTITUENTS)/spec_heat_capacities_p_gas(j - NO_OF_CONDENSED_CONSTITUENTS));
					specific_entropy = spec_heat_capacities_p_gas(j - NO_OF_CONDENSED_CONSTITUENTS)*log(pot_temp);
					init_state -> entropy_densities[j*NO_OF_SCALARS + i] = init_state -> mass_densities[j*NO_OF_SCALARS + i]*specific_entropy;
		    	}
		    }
	    }
	    init_state -> condensed_density_temperatures[i] = solid_water_density[i]*solid_water_temperature[i];
	    init_state -> condensed_density_temperatures[NO_OF_SCALARS + i] = liquid_water_density[i]*liquid_water_temperature[i];
	}
    
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        init_state -> velocity_gas[i] = wind[i];
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
