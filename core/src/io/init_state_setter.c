#include <stdlib.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "io.h"
#include "/home/max/custom_builds/eccodes/include/eccodes.h"
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}

int set_init_data(char FILE_NAME[], State *init_state)
{
    const int START_SECTION = 4;
    long unsigned int no_scalars_h = NUMBER_OF_SCALARS_H;
    long unsigned int no_vectors_h = NUMBER_OF_VECTORS_H;
    FILE *IN_FILE;
    int err = 0;
    codes_handle *handle_pot_temperature = NULL;
    codes_handle *handle_density = NULL;
    codes_handle *handle_wind_h = NULL;
    codes_handle *handle_wind_v = NULL;
    double *pot_temp, *rho, *wind_h, *wind_v;
    pot_temp = malloc(sizeof(double)*NUMBER_OF_SCALARS_H);
    rho = malloc(sizeof(double)*NUMBER_OF_SCALARS_H);
    wind_h = malloc(sizeof(double)*NUMBER_OF_VECTORS_H);
    wind_v = malloc(sizeof(double)*NUMBER_OF_SCALARS_H);
    int retval;
    IN_FILE = fopen(FILE_NAME, "r");
    for (int i = 0; i < NUMBER_OF_LAYERS; i++)
    {
        handle_pot_temperature = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if (retval = codes_get_double_array(handle_pot_temperature, "values", pot_temp, &no_scalars_h))
            ECCERR(retval);
        codes_handle_delete(handle_pot_temperature);
        handle_density = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if (retval = codes_get_double_array(handle_density, "values", rho, &no_scalars_h))
            ECCERR(retval);
        codes_handle_delete(handle_density);
        for (int j = 0; j < NUMBER_OF_SCALARS_H; j++)
        {
            init_state -> density[j + i*NUMBER_OF_SCALARS_H] = rho[j];
            init_state -> density_pot_temp[j + i*NUMBER_OF_SCALARS_H] = rho[j]*pot_temp[j];
        }
        handle_wind_h = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if (retval = codes_get_double_array(handle_wind_h, "values", wind_h, &no_vectors_h))
            ECCERR(retval);
        codes_handle_delete(handle_wind_h);
        for (int j = 0; j < NUMBER_OF_VECTORS_H; j++)
            init_state -> wind[j + i*NUMBER_OF_VECTORS_H + (i + 1)*NUMBER_OF_SCALARS_H] = wind_h[j];
    }
    if (err != 0)
        ECCERR(err);
    for (int i = 0; i < NUMBER_OF_LEVELS; i++)
    {
        handle_wind_v = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if (retval = codes_get_double_array(handle_wind_v, "values", wind_v, &no_scalars_h))
            ECCERR(retval);
        codes_handle_delete(handle_wind_v);
        for (int j = 0; j < NUMBER_OF_SCALARS_H; j++)
            init_state -> wind[j + i*(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H)] = wind_v[j];
    }
    fclose(IN_FILE);
    free(pot_temp);
    free(rho);
    free(wind_h);
    free(wind_v);
    return 0;
}
