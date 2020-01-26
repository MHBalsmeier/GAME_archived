#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <math.h>
#include "/home/max/source/eccodes/include/eccodes.h"
#include "../../models/ess0.0.0/src/enum_and_typedefs.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define FILE_NAME_TRAFO "../../models_geo_prep/calchexgeo0.0.0/nc_files/sphere_res_2_oro_0_rect_trafo.nc"

int main(int argc, char *argv[])
{
   const int productionStatusOfProcessedData = 1;
   const int significanceOfReferenceTime = 1;
   // SETTING THE TYPEOFPROCESSEDDATA
   const int typeOfProcessedData = 1;
   const int START_SECTION = 4;
   char *INPUT_FILE = argv[1];
   char *OUTPUT_FILE = argv[2];
   int STEP;
   sscanf(argv[3], "%d", &STEP);
   int ncid, nc_trafo_id, pressure_id, temperature_id, rho_id, wind_id;
   double *pressure, *temperature, *rho, *wind;
   pressure = (double *) calloc(NUMBER_OF_SCALARS,sizeof(double));
   temperature = (double *) calloc(NUMBER_OF_SCALARS,sizeof(double));
   rho = (double *) calloc(NUMBER_OF_SCALARS,sizeof(double));
   wind = (double *) calloc(NUMBER_OF_VECTORS,sizeof(double));
   int retval;
   if ((retval = nc_open(INPUT_FILE, NC_NOWRITE, &ncid)))
      ERR(retval);
   if ((retval = nc_inq_varid(ncid, "wind", &wind_id)))
      ERR(retval);
   if ((retval = nc_get_var_double(ncid, wind_id, &wind[0])))
      ERR(retval);
   if ((retval = nc_inq_varid(ncid, "pressure", &pressure_id)))
      ERR(retval);
   if ((retval = nc_inq_varid(ncid, "density", &rho_id)))
      ERR(retval);
   if ((retval = nc_inq_varid(ncid, "temperature", &temperature_id)))
      ERR(retval);
   if ((retval = nc_get_var_double(ncid, pressure_id, &pressure[0])))
      ERR(retval);
   if ((retval = nc_get_var_double(ncid, rho_id, &rho[0])))
      ERR(retval);
   if ((retval = nc_get_var_double(ncid, temperature_id, &temperature[0])))
      ERR(retval);
   if ((retval = nc_close(ncid)))
      ERR(retval);
   const int NUMBER_OF_LINES = 0.5*sqrt(NUMBER_OF_SCALARS_H*M_PI);
   const int NUMBER_OF_COLUMNS = 4*NUMBER_OF_LINES/M_PI;
   const double SCALE_HEIGHT = 8e3;
   const double ATMOS_HEIGHT = SCALE_HEIGHT*log(1+NUMBER_OF_LAYERS);
   double lat_deg_inc = NUMBER_OF_LINES/180.0;
   double lon_deg_inc = NUMBER_OF_COLUMNS/360.0;
   double lat_deg_first_grid = -90 + 180.0/(2*NUMBER_OF_LINES);
   double lon_deg_first_grid = lon_deg_inc/2; 
   const int NUMBER_OF_VALUES = NUMBER_OF_LINES*NUMBER_OF_COLUMNS;
   if ((retval = nc_open(FILE_NAME_TRAFO, NC_NOWRITE, &nc_trafo_id)))
      ERR(retval);
   int scalar_hex2rect_indices_id, w_wind_hex2rect_indices_id, wind_hex2rect_indices_id, recov_dimid_id, wind_hex2rect_weights_id, directions_hex2rect_id;
    long *scalar_hex2rect_indices, *w_wind_hex2rect_indices, *wind_hex2rect_indices;
    double *wind_hex2rect_weights, *directions_hex2rect;
    scalar_hex2rect_indices = (long *) calloc(NUMBER_OF_VALUES, sizeof(long));
    w_wind_hex2rect_indices = (long *) calloc(NUMBER_OF_VALUES, sizeof(long));
    wind_hex2rect_indices = (long *) calloc(11*NUMBER_OF_VALUES, sizeof(long));
    wind_hex2rect_weights = (double *) calloc(11*NUMBER_OF_VALUES, sizeof(double));
    directions_hex2rect = (double *) calloc(NUMBER_OF_VALUES, sizeof(double));
    if ((retval = nc_inq_varid(nc_trafo_id, "scalar_indices", &scalar_hex2rect_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(nc_trafo_id, "w_wind_indices", &w_wind_hex2rect_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(nc_trafo_id, "wind_hex2rect_indices", &wind_hex2rect_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(nc_trafo_id, "wind_hex2rect_weights", &wind_hex2rect_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(nc_trafo_id, "directions", &directions_hex2rect_id)))
        ERR(retval);
    if ((retval = nc_get_var_long(nc_trafo_id, scalar_hex2rect_indices_id, &scalar_hex2rect_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(nc_trafo_id, w_wind_hex2rect_indices_id, &w_wind_hex2rect_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_long(nc_trafo_id, wind_hex2rect_indices_id, &wind_hex2rect_indices[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(nc_trafo_id, wind_hex2rect_weights_id, &wind_hex2rect_weights[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(nc_trafo_id, directions_hex2rect_id, &directions_hex2rect[0])))
        ERR(retval);
    if ((retval = nc_close(ncid)))
        ERR(retval);
    double *pressure_rect, *temperature_rect, *rho_rect,*u_wind_rect, *v_wind_rect, *w_wind_rect;
    pressure_rect = (double *) calloc(NUMBER_OF_VALUES,sizeof(double));
    temperature_rect = (double *) calloc(NUMBER_OF_VALUES,sizeof(double));
    rho_rect = (double *) calloc(NUMBER_OF_VALUES,sizeof(double));
    u_wind_rect = (double *) calloc(NUMBER_OF_VALUES,sizeof(double));
    v_wind_rect = (double *) calloc(NUMBER_OF_VALUES,sizeof(double));
    w_wind_rect = (double *) calloc(NUMBER_OF_VALUES,sizeof(double));
    char *SAMPLE_FILE = "test_res_2_sphere_oro_0_scenario_0_template.grb2";
    FILE *SAMPLE;
    SAMPLE = fopen(SAMPLE_FILE,"r");
    int err = 0;
    FILE *OUT_GRIB;
    OUT_GRIB = fopen(OUTPUT_FILE, "w");
    grib_multi_handle *multi_handle=NULL;
    multi_handle=codes_grib_multi_handle_new(NULL);
    double sigma;
    codes_handle *handle_pressure = NULL;
    codes_handle *handle_temperature = NULL;
    codes_handle *handle_rho = NULL;
    codes_handle *handle_u_wind = NULL;
    codes_handle *handle_v_wind = NULL;
    codes_handle *handle_w_wind = NULL;
    for (int layer_index = 0; layer_index < NUMBER_OF_LAYERS; layer_index++)
    {
        sigma = 0.5*(SCALE_HEIGHT/ATMOS_HEIGHT)*(log((1.0+NUMBER_OF_LAYERS)/(layer_index+1))+log((1.0+NUMBER_OF_LAYERS)/(layer_index+2)));
    double ortho_wind, para_wind;
    for (int i = 0; i < NUMBER_OF_LINES; i++)
    {
            for (int j = 0; j < NUMBER_OF_COLUMNS; j++)
            {
                pressure_rect[i*NUMBER_OF_COLUMNS+j] = pressure[scalar_hex2rect_indices[i*NUMBER_OF_COLUMNS+j]];
                temperature_rect[i*NUMBER_OF_COLUMNS+j] = temperature[scalar_hex2rect_indices[i*NUMBER_OF_COLUMNS+j]];
                rho_rect[i*NUMBER_OF_COLUMNS+j] = rho[scalar_hex2rect_indices[i*NUMBER_OF_COLUMNS+j]];
                ortho_wind = wind[wind_hex2rect_indices[i*NUMBER_OF_COLUMNS+j]];
                para_wind = 0;
                for (int k = 0; k < 11; k++)
                {
                    para_wind = para_wind + wind_hex2rect_weights[11*(i*NUMBER_OF_COLUMNS+j)+k]*wind[wind_hex2rect_indices[11*(i*NUMBER_OF_COLUMNS+j)+k]];
                }
                u_wind_rect[i*NUMBER_OF_COLUMNS+j] = ortho_wind*sin(directions_hex2rect[i*NUMBER_OF_COLUMNS+j]) + para_wind*sin(directions_hex2rect[i*NUMBER_OF_COLUMNS+j] + M_PI/2);
                v_wind_rect[i*NUMBER_OF_COLUMNS+j] = ortho_wind*cos(directions_hex2rect[i*NUMBER_OF_COLUMNS+j]) + para_wind*cos(directions_hex2rect[i*NUMBER_OF_COLUMNS+j] + M_PI/2);
                w_wind_rect[i*NUMBER_OF_COLUMNS+j] = 0;
            }
    }
    handle_pressure = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    handle_temperature = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    handle_rho = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    handle_u_wind = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    handle_v_wind = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    handle_w_wind = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    if (retval = codes_set_double_array(handle_pressure,"values",pressure_rect,NUMBER_OF_VALUES))
        printf("Error setting the pressure array to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double_array(handle_temperature,"values",temperature_rect,NUMBER_OF_VALUES))
        printf("Error setting the temperature array to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double_array(handle_rho,"values",rho_rect,NUMBER_OF_VALUES))
        printf("Error setting the density array to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double_array(handle_u_wind,"values",u_wind_rect,NUMBER_OF_VALUES))
        printf("Error setting the u wind array to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double_array(handle_v_wind,"values",v_wind_rect,NUMBER_OF_VALUES))
        printf("Error setting the v wind array to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double_array(handle_w_wind,"values",w_wind_rect,NUMBER_OF_VALUES))
        printf("Error setting the w wind array to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_pressure, "centre", 255))
        printf("Error setting the pressure centre to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "centre", 255))
        printf("Error setting the temperature centre to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "centre", 255))
        printf("Error setting the density centre to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "centre", 255))
        printf("Error setting the u wind centre to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "centre", 255))
        printf("Error setting the v wind centre to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "centre", 255))
        printf("Error setting the w wind centre to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    // setting the productionStatusOfProcessedData to the handles and checking for errors and returning error codes
    if (retval = codes_set_long(handle_pressure, "productionStatusOfProcessedData", productionStatusOfProcessedData))
        printf("Error setting the pressure productionStatusOfProcessedData to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "productionStatusOfProcessedData", productionStatusOfProcessedData))
        printf("Error setting the temperature productionStatusOfProcessedData to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "productionStatusOfProcessedData", productionStatusOfProcessedData))
        printf("Error setting the density productionStatusOfProcessedData to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "productionStatusOfProcessedData", productionStatusOfProcessedData))
        printf("Error setting the u wind productionStatusOfProcessedData to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "productionStatusOfProcessedData", productionStatusOfProcessedData))
        printf("Error setting the v wind productionStatusOfProcessedData to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "productionStatusOfProcessedData", productionStatusOfProcessedData))
        printf("Error setting the w wind productionStatusOfProcessedData to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    // setting the significanceOfReferenceTime to the handles and checking for errors and returning error codes
    if (retval = codes_set_long(handle_pressure, "significanceOfReferenceTime", significanceOfReferenceTime))
        printf("Error setting the pressure significanceOfReferenceTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "significanceOfReferenceTime", significanceOfReferenceTime))
        printf("Error setting the temperature significanceOfReferenceTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "significanceOfReferenceTime", significanceOfReferenceTime))
        printf("Error setting the density significanceOfReferenceTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "significanceOfReferenceTime", significanceOfReferenceTime))
        printf("Error setting the u wind significanceOfReferenceTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "significanceOfReferenceTime", significanceOfReferenceTime))
        printf("Error setting the v wind significanceOfReferenceTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "significanceOfReferenceTime", significanceOfReferenceTime))
        printf("Error setting the w wind significanceOfReferenceTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    // setting the significanceOfReferenceTime to the handles and checking for errors and returning error codes
    if (retval = codes_set_long(handle_pressure, "typeOfProcessedData", typeOfProcessedData))
        printf("Error setting the pressure typeOfProcessedData to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "typeOfProcessedData", typeOfProcessedData))
        printf("Error setting the temperature typeOfProcessedData to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "typeOfProcessedData", typeOfProcessedData))
        printf("Error setting the density typeOfProcessedData to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "typeOfProcessedData", typeOfProcessedData))
        printf("Error setting the u wind typeOfProcessedData to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "typeOfProcessedData", typeOfProcessedData))
        printf("Error setting the v wind typeOfProcessedData to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "typeOfProcessedData", typeOfProcessedData))
        printf("Error setting the w wind typeOfProcessedData to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    // setting the forecastTime to the handles and checking for errors and returning error codes
    if (retval = codes_set_long(handle_pressure, "forecastTime", STEP))
        printf("Error setting the pressure forecastTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "forecastTime", STEP))
        printf("Error setting the temperature forecastTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "forecastTime", STEP))
        printf("Error setting the density forecastTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "forecastTime", STEP))
        printf("Error setting the u wind forecastTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "forecastTime", STEP))
        printf("Error setting the v wind forecastTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "forecastTime", STEP))
        printf("Error setting the w wind forecastTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    // setting the typeOfFirstFixedSurface to the handles and checking for errors and returning error codes
    if (retval = codes_set_long(handle_pressure, "typeOfFirstFixedSurface", 104))
        printf("Error setting the pressure typeOfFirstFixedSurface to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "typeOfFirstFixedSurface", 104))
        printf("Error setting the temperature typeOfFirstFixedSurface to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "typeOfFirstFixedSurface", 104))
        printf("Error setting the density typeOfFirstFixedSurface to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "typeOfFirstFixedSurface", 104))
        printf("Error setting the u wind typeOfFirstFixedSurface to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "typeOfFirstFixedSurface", 104))
        printf("Error setting the v wind typeOfFirstFixedSurface to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "typeOfFirstFixedSurface", 104))
        printf("Error setting the w wind typeOfFirstFixedSurface to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_pressure, "level", sigma*10e5))
        printf("Error setting the pressure level to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "level", sigma*10e5))
        printf("Error setting the temperature level to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "level", sigma*10e5))
        printf("Error setting the density level to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "level", sigma*10e5))
        printf("Error setting the u wind level to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "level", sigma*10e5))
        printf("Error setting the v wind level to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "level", sigma*10e5))
        printf("Error setting the w wind level to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_pressure, "jScansPositively", 1))
        printf("Error setting the pressure jScansPositively to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "jScansPositively", 1))
        printf("Error setting the temperature jScansPositively to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "jScansPositively", 1))
        printf("Error setting the density jScansPositively to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "jScansPositively", 1))
        printf("Error setting the u wind jScansPositively to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "jScansPositively", 1))
        printf("Error setting the v wind jScansPositively to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "jScansPositively", 1))
        printf("Error setting the w wind jScansPositively to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_pressure, "jPointsAreConsecutive", 1))
        printf("Error setting the pressure jPointsAreConsecutive to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "jPointsAreConsecutive", 1))
        printf("Error setting the temperature jPointsAreConsecutive to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "jPointsAreConsecutive", 1))
        printf("Error setting the density jPointsAreConsecutive to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "jPointsAreConsecutive", 1))
        printf("Error setting the u wind jPointsAreConsecutive to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "jPointsAreConsecutive", 1))
        printf("Error setting the v wind jPointsAreConsecutive to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "jPointsAreConsecutive", 1))
        printf("Error setting the w wind jPointsAreConsecutive to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    // setting the jDirectionIncrementInDegrees to the handles and checking for errors and returning error codes
    if (retval = codes_set_double(handle_pressure, "jDirectionIncrementInDegrees", lat_deg_inc))
        printf("Error setting the pressure jDirectionIncrementInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_temperature, "jDirectionIncrementInDegrees", lat_deg_inc))
        printf("Error setting the temperature jDirectionIncrementInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_rho, "jDirectionIncrementInDegrees", lat_deg_inc))
        printf("Error setting the density jDirectionIncrementInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_u_wind, "jDirectionIncrementInDegrees", lat_deg_inc))
        printf("Error setting the u wind jDirectionIncrementInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_v_wind, "jDirectionIncrementInDegrees", lat_deg_inc))
        printf("Error setting the v wind jDirectionIncrementInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_w_wind, "jDirectionIncrementInDegrees", lat_deg_inc))
        printf("Error setting the w wind jDirectionIncrementInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_pressure, "iDirectionIncrementInDegrees", lon_deg_inc))
        printf("Error setting the pressure iDirectionIncrementInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_temperature, "iDirectionIncrementInDegrees", lon_deg_inc))
        printf("Error setting the temperature iDirectionIncrementInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_rho, "iDirectionIncrementInDegrees", lon_deg_inc))
        printf("Error setting the density iDirectionIncrementInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_u_wind, "iDirectionIncrementInDegrees", lon_deg_inc))
        printf("Error setting the u wind iDirectionIncrementInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_v_wind, "iDirectionIncrementInDegrees", lon_deg_inc))
        printf("Error setting the v wind iDirectionIncrementInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_w_wind, "iDirectionIncrementInDegrees", lon_deg_inc))
        printf("Error setting the w wind iDirectionIncrementInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_pressure, "latitudeOfFirstGridPointInDegrees", lat_deg_first_grid))
        printf("Error setting the pressure latitudeOfFirstGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_temperature, "latitudeOfFirstGridPointInDegrees", lat_deg_first_grid))
        printf("Error setting the temperature latitudeOfFirstGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_rho, "latitudeOfFirstGridPointInDegrees", lat_deg_first_grid))
        printf("Error setting the density latitudeOfFirstGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_u_wind, "latitudeOfFirstGridPointInDegrees", lat_deg_first_grid))
        printf("Error setting the u wind latitudeOfFirstGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_v_wind, "latitudeOfFirstGridPointInDegrees", lat_deg_first_grid))
        printf("Error setting the v wind latitudeOfFirstGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_w_wind, "latitudeOfFirstGridPointInDegrees", lat_deg_first_grid))
        printf("Error setting the w wind latitudeOfFirstGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_pressure, "latitudeOfLastGridPointInDegrees", -lat_deg_first_grid))
        printf("Error setting the pressure latitudeOfLastGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_temperature, "latitudeOfLastGridPointInDegrees", -lat_deg_first_grid))
        printf("Error setting the temperature latitudeOfLastGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_rho, "latitudeOfLastGridPointInDegrees", -lat_deg_first_grid))
        printf("Error setting the density latitudeOfLastGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_u_wind, "latitudeOfLastGridPointInDegrees", -lat_deg_first_grid))
        printf("Error setting the u wind latitudeOfLastGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_v_wind, "latitudeOfLastGridPointInDegrees", -lat_deg_first_grid))
        printf("Error setting the v wind latitudeOfLastGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_w_wind, "latitudeOfLastGridPointInDegrees", -lat_deg_first_grid))
        printf("Error setting the w wind latitudeOfLastGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_pressure, "longitudeOfFirstGridPointInDegrees", lon_deg_first_grid))
        printf("Error setting the pressure longitudeOfFirstGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_temperature, "longitudeOfFirstGridPointInDegrees", lon_deg_first_grid))
        printf("Error setting the temperature longitudeOfFirstGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_rho, "longitudeOfFirstGridPointInDegrees", lon_deg_first_grid))
        printf("Error setting the density longitudeOfFirstGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_u_wind, "longitudeOfFirstGridPointInDegrees", lon_deg_first_grid))
        printf("Error setting the u wind longitudeOfFirstGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_v_wind, "longitudeOfFirstGridPointInDegrees", lon_deg_first_grid))
        printf("Error setting the v wind longitudeOfFirstGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_w_wind, "longitudeOfFirstGridPointInDegrees", lon_deg_first_grid))
        printf("Error setting the w wind longitudeOfFirstGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_pressure, "longitudeOfLastGridPointInDegrees", 360-lon_deg_first_grid))
        printf("Error setting the pressure longitudeOfLastGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_temperature, "longitudeOfLastGridPointInDegrees", 360-lon_deg_first_grid))
        printf("Error setting the temperature longitudeOfLastGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_rho, "longitudeOfLastGridPointInDegrees", 360-lon_deg_first_grid))
        printf("Error setting the density longitudeOfLastGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_u_wind, "longitudeOfLastGridPointInDegrees", 360-lon_deg_first_grid))
        printf("Error setting the u wind longitudeOfLastGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_v_wind, "longitudeOfLastGridPointInDegrees", 360-lon_deg_first_grid))
        printf("Error setting the v wind longitudeOfLastGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_double(handle_w_wind, "longitudeOfLastGridPointInDegrees", 360-lon_deg_first_grid))
        printf("Error setting the w wind longitudeOfLastGridPointInDegrees to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_pressure, "Ni", NUMBER_OF_COLUMNS))
        printf("Error setting the pressure wind Ni to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "Ni", NUMBER_OF_COLUMNS))
        printf("Error setting the temperature wind Ni to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "Ni", NUMBER_OF_COLUMNS))
        printf("Error setting the density wind Ni to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "Ni", NUMBER_OF_COLUMNS))
        printf("Error setting the u wind Ni to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "Ni", NUMBER_OF_COLUMNS))
        printf("Error setting the v wind Ni to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "Ni", NUMBER_OF_COLUMNS))
        printf("Error setting the w wind Ni to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_pressure, "Nj", NUMBER_OF_LINES))
        printf("Error setting the pressure Nj to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "Nj", NUMBER_OF_LINES))
        printf("Error setting the temperature Nj to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "Nj", NUMBER_OF_LINES))
        printf("Error setting the density Nj to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "Nj", NUMBER_OF_LINES))
        printf("Error setting the u wind Nj to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "Nj", NUMBER_OF_LINES))
        printf("Error setting the v wind Nj to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "Nj", NUMBER_OF_LINES))
        printf("Error setting the w wind Nj to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_pressure, "dataDate", 20000101))
        printf("Error setting the pressure dataDate to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "dataDate", 20000101))
        printf("Error setting the temperature dataDate to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "dataDate", 20000101))
        printf("Error setting the density dataDate to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "dataDate", 20000101))
        printf("Error setting the u wind dataDate to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "dataDate", 20000101))
        printf("Error setting the v wind dataDate to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "dataDate", 20000101))
        printf("Error setting the w wind dataDate to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_pressure, "dataTime", 0000))
        printf("Error setting the pressure dataTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "dataTime", 0000))
        printf("Error setting the temperature dataTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "dataTime", 0000))
        printf("Error setting the density dataTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "dataTime", 0000))
        printf("Error setting the u wind dataTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "dataTime", 0000))
        printf("Error setting the v wind dataTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "dataTime", 0000))
        printf("Error setting the w wind dataTime to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_pressure, "parameterCategory", 3))
        printf("Error setting the pressure parameterCategory to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "parameterCategory", 0))
        printf("Error setting the temperature parameterCategory to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "parameterCategory", 3))
        printf("Error setting the density parameterCategory to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "parameterCategory", 2))
        printf("Error setting the u wind parameterCategory to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "parameterCategory", 2))
        printf("Error setting the v wind parameterCategory to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "parameterCategory", 2))
        printf("Error setting the w wind parameterCategory to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_pressure, "parameterNumber", 0))
        printf("Error setting the pressure parameterNumber to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_temperature, "parameterNumber", 0))
        printf("Error setting the temperature wind parameterNumber to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_rho, "parameterNumber", 10))
        printf("Error setting the density wind parameterNumber to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_u_wind, "parameterNumber", 2))
        printf("Error setting the u wind parameterNumber to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_v_wind, "parameterNumber", 3))
        printf("Error setting the v wind parameterNumber to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    if (retval = codes_set_long(handle_w_wind, "parameterNumber", 9))
        printf("Error setting the w wind parameterNumber to the respective handle in File %s. Error code: %d.\n",OUTPUT_FILE,retval);
    codes_grib_multi_handle_append(handle_pressure, START_SECTION, multi_handle);
    codes_grib_multi_handle_append(handle_temperature, START_SECTION, multi_handle);
    codes_grib_multi_handle_append(handle_rho, START_SECTION, multi_handle);
    codes_grib_multi_handle_append(handle_u_wind, START_SECTION, multi_handle);
    codes_grib_multi_handle_append(handle_v_wind, START_SECTION, multi_handle);
    codes_grib_multi_handle_append(handle_w_wind, START_SECTION, multi_handle);
   }
   codes_handle_delete(handle_pressure);
   codes_handle_delete(handle_temperature);
   codes_handle_delete(handle_rho);
   codes_handle_delete(handle_u_wind);
   codes_handle_delete(handle_v_wind);
   codes_handle_delete(handle_w_wind);
   codes_grib_multi_handle_write(multi_handle, OUT_GRIB);
   codes_grib_multi_handle_delete(multi_handle);
   fclose(SAMPLE);
   fclose(OUT_GRIB);
   free(pressure_rect);
   free(temperature_rect);
   free(rho_rect);
   free(u_wind_rect);
   free(v_wind_rect);
   free(w_wind_rect);
   free(wind);
   free(rho);
   free(temperature);
   free(pressure);
   return 0;
}
