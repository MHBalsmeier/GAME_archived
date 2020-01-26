#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "../enum_and_typedefs.h"
#include "io.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

void set_init_data(char FILE_NAME[], State init_state)
{
    int ncid, temperature_id, rho_id, wind_id;
    double temperature[NUMBER_OF_SCALARS];
    double rho[NUMBER_OF_SCALARS];
    double wind[NUMBER_OF_VECTORS];
    int retval;
    double pressure;
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature", &temperature_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "density", &rho_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "wind", &wind_id)))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, temperature_id, &temperature[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, rho_id, &rho[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, wind_id, &wind[0])))
        ERR(retval);
    for (int i = 0; i < NUMBER_OF_SCALARS; i++)
    {
        init_state.density[i] = rho[i];
        init_state.pressure[i] = rho[i]*R_D*temperature[i];
        init_state.pot_temp[i] = temperature[i]*pow(P_0/init_state.pressure[i], R_D/C_P);
    }
    for (int i = 0; i < NUMBER_OF_VECTORS; i++)
        init_state.wind[i] = wind[i];
    if ((retval = nc_close(ncid)))
        ERR(retval);
}
