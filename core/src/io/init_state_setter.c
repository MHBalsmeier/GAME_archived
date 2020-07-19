/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
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
    double *pot_temp = malloc(NO_OF_SCALARS*sizeof(double));
    double *rho = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *water_vapour_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_temperature = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_temperature = malloc(NO_OF_SCALARS*sizeof(double));
    int retval, ncid;
    double R_h, c_h_v, density_h_micro;
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        NCERR(retval);
    int pot_temp_id, density_dry_id, wind_id, density_vapour_id, density_liquid_id, density_solid_id, temperature_liquid_id, temperature_solid_id;
    if ((retval = nc_inq_varid(ncid, "potential_temperature",&pot_temp_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_dry", &density_dry_id)))
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
    if ((retval = nc_get_var_double(ncid, pot_temp_id, &pot_temp[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, density_dry_id, &rho[0])))
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
        init_state -> density_dry[i] = rho[i];
        R_h = gas_constant_diagnostics(rho[i], water_vapour_density[i]);
        c_h_v = spec_heat_cap_diagnostics_p(rho[i], water_vapour_density[i]);
        density_h_micro = calc_micro_density(rho[i] + water_vapour_density[i], solid_water_density[i] + liquid_water_density[i]);
        init_state -> t_tilde[i] = init_state -> density_dry[i]*pow(pot_temp[i], C_D_P/C_D_V)*pow(R_D*init_state -> density_dry[i]/P_0, R_D/C_D_V);
        init_state -> entropy_gas[i] = rho[i]*(C_D_P*log(pot_temp[i]) + entropy_constant_d) + water_vapour_density[i]*(C_V_P*log(pot_temp[i]) + M_D/M_V*DELTA_C_V_P*R_h/c_h_v*log(R_h*pot_temp[i]*density_h_micro/P_0) + entropy_constant_d);
        if (NO_OF_TRACERS > 0)
        {
            init_state -> tracer_densities[i] = solid_water_density[i];
            init_state -> tracer_densities[NO_OF_SCALARS + i] = liquid_water_density[i];
            init_state -> tracer_densities[2*NO_OF_SCALARS + i] = water_vapour_density[i];
            init_state -> tracer_density_temperatures[i] = solid_water_density[i]*solid_water_temperature[i];
            init_state -> tracer_density_temperatures[NO_OF_SCALARS + i] = liquid_water_density[i]*liquid_water_temperature[i];
        }
    }
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        init_state -> velocity_gas[i] = wind[i];
    }
    free(pot_temp);
    free(rho);
    free(wind);
    free(water_vapour_density);
    free(liquid_water_density);
    free(solid_water_density);
    free(liquid_water_temperature);
    free(solid_water_temperature);
    return 0;
}
