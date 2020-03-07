#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include "../enum_and_typedefs.h"
#include "io.h"
#include "../diagnostics/diagnostics.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

void write_out(State state_write_out, double t_init, double t_write, int write_out_index, char output_foldername[])
{
    char add_time_str[write_out_index/10+1];
    char pre_string[] = "init+";
    char append_string[] = "h.nc";
    sprintf(add_time_str, "%d", write_out_index);
    int out_file_length = strlen(output_foldername) + 1 + strlen(pre_string) + write_out_index/10 + 1 + strlen(append_string);
    char FILE_WRITE[out_file_length];
    for (int i = 0; i < out_file_length + 1; i++)
    {
        if (i < strlen(output_foldername))
            FILE_WRITE[i] = output_foldername[i];
        if (i == strlen(output_foldername))
            FILE_WRITE[i] = '/';
        if (i >= strlen(output_foldername) + 1 && i < strlen(output_foldername) + 1 + strlen(pre_string))
            FILE_WRITE[i] = pre_string[i - (strlen(output_foldername) + 1)];
        if (i >= strlen(output_foldername) + 1 + strlen(pre_string) && i < out_file_length - strlen(append_string))
            FILE_WRITE[i] = add_time_str[i - (strlen(output_foldername) + 1 + strlen(pre_string))];
        if (i >= out_file_length - strlen(append_string))
            FILE_WRITE[i] = append_string[i - (out_file_length - strlen(append_string))];
    }
    Scalar_field pressure, temperature;
    pressure_diagnostics(state_write_out.pot_temp, state_write_out.density, pressure);
    temperature_diagnostics(state_write_out.pot_temp, pressure, temperature);
    int ncid, scalar_dimid, vector_dimid, var_dimid, pressure_id, rho_id, temperature_id, wind_id;
    int dimids_scalar[2];
    int scalar_index,  retval;
    if ((retval = nc_create(FILE_WRITE, NC_CLOBBER, &ncid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "scalar_index", NUMBER_OF_SCALARS, &scalar_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "vector_index", NUMBER_OF_VECTORS, &vector_dimid)))
        ERR(retval);
    if ((retval = nc_def_var(ncid, "pressure", NC_DOUBLE, 1, &scalar_dimid, &pressure_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid, pressure_id, "units", strlen("hPa"), "hPa")))
        ERR(retval);
    if ((retval = nc_def_var(ncid, "density", NC_DOUBLE, 1, &scalar_dimid, &rho_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid, rho_id, "units", strlen("kg/m^3"), "kg/m^3")))
        ERR(retval);
    if ((retval = nc_def_var(ncid, "temperature", NC_DOUBLE, 1, &scalar_dimid, &temperature_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid, temperature_id, "units", strlen("K"), "K")))
        ERR(retval);
    if ((retval = nc_def_var(ncid, "wind", NC_DOUBLE, 1, &vector_dimid, &wind_id)))
        ERR(retval);
    if ((retval = nc_put_att_text(ncid, wind_id, "units", strlen("m/s"), "m/s")))
        ERR(retval);
    if ((retval = nc_enddef(ncid)))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid, pressure_id, &pressure[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid, rho_id, &state_write_out.density[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid, temperature_id, &temperature[0])))
        ERR(retval);
    if ((retval = nc_put_var_double(ncid, wind_id, &state_write_out.wind[0])))
        ERR(retval);
    if ((retval = nc_close(ncid)))
        ERR(retval);
}
