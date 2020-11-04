/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include <stdio.h>
#include <stdlib.h>
#include "../enum_and_typedefs.h"
#include "io.h"
#include "../diagnostics/diagnostics.h"
#include <netcdf.h>
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

int set_init_data(char FILE_NAME[], State *init_state)
{
    double *density_dry = malloc(NO_OF_SCALARS*sizeof(double));
    double *temperature_gas = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *pot_temperature = malloc(NO_OF_SCALARS*sizeof(double));
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
        init_state -> density_dry[i] = density_dry[i];
        init_state -> temperature_gas[i] = temperature_gas[i];
    }
    if (NO_OF_TRACERS > 0)
    {
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
	        init_state -> tracer_densities[i] = solid_water_density[i];
	        init_state -> tracer_densities[NO_OF_SCALARS + i] = liquid_water_density[i];
	        init_state -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i] = water_vapour_density[i];
		}
    }    
	pot_temp_diagnostics(init_state, pot_temperature);
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        init_state -> entropy_density_dry[i] = density_dry[i]*(C_D_P*log(pot_temperature[i]) + ENTROPY_CONSTANT_D);
    }
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        init_state -> velocity_gas[i] = wind[i];
    }
    if (NO_OF_TRACERS > 0)
    {
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
	        init_state -> tracer_entropy_densities[i] = init_state -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i]*(C_V_P*log(pot_temperature[i]) + ENTROPY_CONSTANT_V);
	        init_state -> tracer_density_temperatures[i] = solid_water_density[i]*solid_water_temperature[i];
	        init_state -> tracer_density_temperatures[NO_OF_SCALARS + i] = liquid_water_density[i]*liquid_water_temperature[i];
		}
    }
    free(pot_temperature);
    free(density_dry);
    free(wind);
    free(water_vapour_density);
    free(liquid_water_density);
    free(solid_water_density);
    free(liquid_water_temperature);
    free(solid_water_temperature);
    return 0;
}
