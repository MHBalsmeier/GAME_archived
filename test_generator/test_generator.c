#include "../models/ess0.0.0/src/enum_and_typedefs.h"
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#define FILE_NAME "nc_files/test_res_2_sphere_oro_0_scenario_0.nc"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int main(int argc, char *argv[])
{
   int ncid, scalar_dimid, vector_dimid, var_dimid, pressure_id, rho_id, temperature_id, wind_id;
   double *pressure, *temperature, *rho, *wind;
   pressure = malloc(sizeof(double)*NUMBER_OF_SCALARS);
   temperature = malloc(sizeof(double)*NUMBER_OF_SCALARS);
   rho = malloc(sizeof(double)*NUMBER_OF_SCALARS);
   wind = malloc(sizeof(double)*NUMBER_OF_VECTORS);
   const double G = 9.8;
   const double SCALE_HEIGHT = 8e3;
   const double TROPO_HEIGHT = 12e3;
   const double T_SFC = 273.15+15;
   const double TEMP_GRADIENT = -0.65/100;
   const double TROPO_TEMP = T_SFC+TROPO_HEIGHT*TEMP_GRADIENT;
   const double ATMOS_HEIGHT = SCALE_HEIGHT*log(1+NUMBER_OF_LAYERS);
   int scalar_index,  retval, layer_index;
   double sigma, z_height;
   for (int i = 0; i < NUMBER_OF_SCALARS; i++)
   {
      layer_index = i/NUMBER_OF_SCALARS_H;
      sigma = 0.5*(SCALE_HEIGHT/ATMOS_HEIGHT)*(log((1.0+NUMBER_OF_LAYERS)/(layer_index+1))+log((1.0+NUMBER_OF_LAYERS)/(layer_index+2)));
      z_height = ATMOS_HEIGHT*sigma;
      if (z_height < TROPO_HEIGHT)
      {
        temperature[i] = T_SFC+z_height*TEMP_GRADIENT;
        pressure[i] = P_0*pow(1 + TEMP_GRADIENT*z_height/T_SFC,-G/(R_D*TEMP_GRADIENT));
      }
      else
      {
        temperature[i] = TROPO_TEMP;
        pressure[i] = P_0*pow(1 + TEMP_GRADIENT*TROPO_HEIGHT/T_SFC,-G/(R_D*TEMP_GRADIENT))*exp(-G*(z_height-TROPO_HEIGHT)/(R_D*TROPO_TEMP));
      }
      rho[i] = pressure[i]/(R_D*temperature[i]);
   }
   for (int i = 0; i < NUMBER_OF_VECTORS; i++)
   {
       wind[i] = 0;
   }
    if ((retval = nc_create(FILE_NAME, NC_CLOBBER, &ncid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "scalar_index", NUMBER_OF_SCALARS, &scalar_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "vector_index", NUMBER_OF_VECTORS, &vector_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "variable_index", 3, &var_dimid)))
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
   if ((retval = nc_put_var_double(ncid, rho_id, &rho[0])))
      ERR(retval);
   if ((retval = nc_put_var_double(ncid, temperature_id, &temperature[0])))
      ERR(retval);
   if ((retval = nc_put_var_double(ncid, wind_id, &wind[0])))
      ERR(retval);
   if ((retval = nc_close(ncid)))
      ERR(retval);
   free(pressure);
   free(temperature);
   free(rho);
   return 0;
}
