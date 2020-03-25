#include "../core/src/enum_and_typedefs.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "/usr/src/eccodes/include/eccodes.h"
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}

int main(int argc, char *argv[])
{
    const int START_SECTION = 4;
    char *SAMPLE_FILE_SCALAR = "grib_files/scalar_field_blueprint_res_id_2.grb2";
    char *SAMPLE_FILE_VECTOR = "grib_files/vector_field_blueprint_res_id_2.grb2";
    FILE *SAMPLE_SCALAR;
    FILE *SAMPLE_VECTOR;
    int err = 0;
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    char *OUTPUT_FILE = "grib_files/test_1_res_4_oro_0.grb2";
    codes_handle *handle_pot_temperature = NULL;
    codes_handle *handle_density = NULL;
    codes_handle *handle_wind_h = NULL;
    codes_handle *handle_wind_v = NULL;
    double *pressure, *pot_temperature, *temperature, *rho, *wind_h, *wind_v;
    pressure = malloc(sizeof(double)*NUMBER_OF_SCALARS_H);
    pot_temperature = malloc(sizeof(double)*NUMBER_OF_SCALARS_H);
    temperature = malloc(sizeof(double)*NUMBER_OF_SCALARS_H);
    rho = malloc(sizeof(double)*NUMBER_OF_SCALARS_H);
    wind_h = malloc(sizeof(double)*NUMBER_OF_VECTORS_H);
    wind_v = malloc(sizeof(double)*NUMBER_OF_SCALARS_H);
    const double G = 9.8;
    const double TROPO_HEIGHT = 12e3;
    const double T_SFC = 273.15 + 15;
    const double TEMP_GRADIENT = -0.65/100;
    const double TROPO_TEMP = T_SFC + TROPO_HEIGHT*TEMP_GRADIENT;
    const double ATMOS_HEIGHT = SCALE_HEIGHT*log(1 + NUMBER_OF_LAYERS);
    int scalar_index,  retval;
    double sigma, z_height;
    handle_pot_temperature = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    handle_density = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    SAMPLE_VECTOR = fopen(SAMPLE_FILE_VECTOR, "r");
    handle_wind_h = codes_handle_new_from_file(NULL, SAMPLE_VECTOR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_VECTOR);
    for (int i = 0; i < NUMBER_OF_LAYERS; i++)
    {
        sigma = 0.5*(SCALE_HEIGHT/ATMOS_HEIGHT)*(log((1.0 + NUMBER_OF_LAYERS)/(i + 1)) + log((1.0 + NUMBER_OF_LAYERS)/(i + 2)));
        z_height = ATMOS_HEIGHT*sigma;
        for (int j = 0; j < NUMBER_OF_SCALARS_H; j++)
        {
            if (z_height < TROPO_HEIGHT)
            {
                temperature[j] = T_SFC + z_height*TEMP_GRADIENT;
                pressure[j] = P_0*pow(1 + TEMP_GRADIENT*z_height/T_SFC, -G/(R_D*TEMP_GRADIENT));
            }
            else
            {
                temperature[j] = TROPO_TEMP;
                pressure[j] = P_0*pow(1 + TEMP_GRADIENT*TROPO_HEIGHT/T_SFC, -G/(R_D*TEMP_GRADIENT))*exp(-G*(z_height - TROPO_HEIGHT)/(R_D*TROPO_TEMP));
            }
            rho[j] = pressure[j]/(R_D*temperature[j]);
            pot_temperature[j] = temperature[j]*pow(P_0/pressure[j], R_D/C_P);
        }
        if (retval = codes_set_long(handle_pot_temperature, "discipline", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "centre", 255))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "significanceOfReferenceTime", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "productionStatusOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "typeOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "indicatorOfUnitOfTimeRange", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "dataDate", 20000101))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "dataTime", 0000))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "typeOfGeneratingProcess", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "parameterCategory", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "parameterNumber", 2))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "typeOfFirstFixedSurface", 102))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "scaleFactorOfFirstFixedSurface", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "scaledValueOfFirstFixedSurface", z_height))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature, "level", i))
            ECCERR(retval);
        if (retval = codes_set_double_array(handle_pot_temperature, "values", pot_temperature, NUMBER_OF_SCALARS_H))
            ECCERR(retval);
        if (i == 0)
            codes_write_message(handle_pot_temperature, OUTPUT_FILE, "w");
        else
            codes_write_message(handle_pot_temperature, OUTPUT_FILE, "a");
        if (retval = codes_set_long(handle_density, "discipline", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "centre", 255))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "significanceOfReferenceTime", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "productionStatusOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "typeOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "indicatorOfUnitOfTimeRange", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "dataDate", 20000101))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "dataTime", 0000))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "typeOfGeneratingProcess", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "parameterCategory", 3))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "parameterNumber", 10))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "typeOfFirstFixedSurface", 102))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "scaleFactorOfFirstFixedSurface", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "scaledValueOfFirstFixedSurface", z_height))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density, "level", i))
            ECCERR(retval);
        if (retval = codes_set_double_array(handle_density, "values", rho, NUMBER_OF_SCALARS_H))
            ECCERR(retval);
        codes_write_message(handle_density, OUTPUT_FILE, "a");
        for (int j = 0; j < NUMBER_OF_VECTORS_H; j++)
            wind_h[j] = 0;
        if (i == NUMBER_OF_LAYERS - 2)
            wind_h[NUMBER_OF_SCALARS_H/2] = 10;
        if (retval = codes_set_long(handle_wind_h, "discipline", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "centre", 255))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "significanceOfReferenceTime", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "productionStatusOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "typeOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "indicatorOfUnitOfTimeRange", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "dataDate", 20000101))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "dataTime", 0000))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "typeOfGeneratingProcess", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "parameterCategory", 2))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "parameterNumber", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "typeOfFirstFixedSurface", 102))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "scaleFactorOfFirstFixedSurface", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "scaledValueOfFirstFixedSurface", z_height))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h, "level", i))
            ECCERR(retval);
        if (retval = codes_set_double_array(handle_wind_h, "values", wind_h, NUMBER_OF_VECTORS_H))
            ECCERR(retval);
        codes_write_message(handle_wind_h, OUTPUT_FILE, "a");
    }
    codes_handle_delete(handle_pot_temperature);
    codes_handle_delete(handle_density);
    codes_handle_delete(handle_wind_h);
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    handle_wind_v = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    for (int i = 0; i < NUMBER_OF_LEVELS; i++)
    {
        sigma = (SCALE_HEIGHT/ATMOS_HEIGHT)*log((1.0 + NUMBER_OF_LAYERS)/(i + 1));
        z_height = ATMOS_HEIGHT*sigma;
        if (retval = codes_set_long(handle_wind_v, "discipline", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "centre", 255))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "significanceOfReferenceTime", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "productionStatusOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "typeOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "indicatorOfUnitOfTimeRange", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "dataDate", 20000101))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "dataTime", 0000))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "typeOfGeneratingProcess", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "parameterCategory", 2))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "parameterNumber", 9))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "typeOfFirstFixedSurface", 102))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "scaleFactorOfFirstFixedSurface", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "scaledValueOfFirstFixedSurface", z_height))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v, "level", i))
            ECCERR(retval);
        if (retval = codes_set_double_array(handle_wind_v, "values", wind_v, NUMBER_OF_SCALARS_H))
            ECCERR(retval);
        codes_write_message(handle_wind_v, OUTPUT_FILE, "a");
    }
    free(pressure);
    free(pot_temperature);
    free(temperature);
    free(rho);
    free(wind_h);
    free(wind_v);
    return 0;
}
