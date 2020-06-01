#include <stdlib.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "io.h"
#include "../diagnostics/diagnostics.h"
#include "eccodes.h"
#include "time00.h"
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}

int set_init_data(char FILE_NAME[], State *init_state, double *t_init, int add_comps_bool)
{
    long unsigned int no_scalars_h = NUMBER_OF_SCALARS_H;
    long unsigned int no_vectors_h = NUMBER_OF_VECTORS_H;
    FILE *IN_FILE;
    int err = 0;
    *t_init = 0;
    codes_handle *handle_pot_temperature = NULL;
    codes_handle *handle_density = NULL;
    codes_handle *handle_wind_h = NULL;
    codes_handle *handle_wind_v = NULL;
    codes_handle *handle_water_vapour_density = NULL;
    codes_handle *handle_liquid_water_density = NULL;
    codes_handle *handle_solid_water_density = NULL;
    codes_handle *handle_liquid_water_temp = NULL;
    codes_handle *handle_solid_water_temp = NULL;
    double *pot_temp = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *rho = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *water_vapour_density = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *liquid_water_density = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *solid_water_density = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *liquid_water_temp = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *solid_water_temp = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *wind_h = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *wind_v = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    int retval;
    long data_date, data_time;
    double R_h, c_h_v, density_h_micro;
    IN_FILE = fopen(FILE_NAME, "r");
    for (int i = 0; i < NUMBER_OF_LAYERS; ++i)
    {
        handle_pot_temperature = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if ((retval = codes_get_double_array(handle_pot_temperature, "values", pot_temp, &no_scalars_h)))
            ECCERR(retval);
        codes_handle_delete(handle_pot_temperature);
        handle_density = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if ((retval = codes_get_double_array(handle_density, "values", rho, &no_scalars_h)))
            ECCERR(retval);
        codes_handle_delete(handle_density);
        handle_water_vapour_density = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if ((retval = codes_get_double_array(handle_water_vapour_density, "values", water_vapour_density, &no_scalars_h)))
            ECCERR(retval);
        codes_handle_delete(handle_water_vapour_density);
        handle_liquid_water_density = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if ((retval = codes_get_double_array(handle_liquid_water_density, "values", liquid_water_density, &no_scalars_h)))
            ECCERR(retval);
        codes_handle_delete(handle_liquid_water_density);
        handle_solid_water_density = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if ((retval = codes_get_double_array(handle_solid_water_density, "values", solid_water_density, &no_scalars_h)))
            ECCERR(retval);
        codes_handle_delete(handle_solid_water_density);
        handle_liquid_water_temp = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if ((retval = codes_get_double_array(handle_liquid_water_temp, "values", liquid_water_temp, &no_scalars_h)))
            ECCERR(retval);
        codes_handle_delete(handle_liquid_water_temp);
        handle_solid_water_temp = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if ((retval = codes_get_double_array(handle_solid_water_temp, "values", solid_water_temp, &no_scalars_h)))
            ECCERR(retval);
        codes_handle_delete(handle_solid_water_temp);
        for (int j = 0; j < NUMBER_OF_SCALARS_H; ++j)
        {
            init_state -> density[j + i*NUMBER_OF_SCALARS_H] = rho[j];
            R_h = gas_constant_diagnostics(rho[j], water_vapour_density[j]);
            c_h_v = spec_heat_cap_diagnostics_p(rho[j], water_vapour_density[j]);
            density_h_micro = calc_micro_density(rho[j] + water_vapour_density[j], solid_water_density[j] + liquid_water_density[j]);
            init_state -> density_entropy[j + i*NUMBER_OF_SCALARS_H] = rho[j]*(C_D_P*log(pot_temp[j]) + entropy_constant_d) + water_vapour_density[j]*(C_V_P*log(pot_temp[j]) + M_D/M_V*DELTA_C_V_P*R_h/c_h_v*log(R_h*pot_temp[j]*density_h_micro/P_0) + entropy_constant_d);
            if (NUMBER_OF_ADD_COMPS > 0)
            {
                init_state -> add_comp_densities[j + i*NUMBER_OF_SCALARS_H] = solid_water_density[j];
                init_state -> add_comp_densities[NUMBER_OF_SCALARS + j + i*NUMBER_OF_SCALARS_H] = liquid_water_density[j];
                init_state -> add_comp_densities[2*NUMBER_OF_SCALARS + j + i*NUMBER_OF_SCALARS_H] = water_vapour_density[j];
                init_state -> add_comp_temps[j + i*NUMBER_OF_SCALARS_H] = solid_water_temp[j];
                init_state -> add_comp_temps[NUMBER_OF_SCALARS + j + i*NUMBER_OF_SCALARS_H] = liquid_water_temp[j];
            }
        }
        handle_wind_h = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if ((retval = codes_get_double_array(handle_wind_h, "values", wind_h, &no_vectors_h)))
            ECCERR(retval);
        for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
            init_state -> wind[j + i*NUMBER_OF_VECTORS_H + (i + 1)*NUMBER_OF_VECTORS_V] = wind_h[j];
        codes_get_long(handle_wind_h, "dataDate", &data_date);
        codes_get_long(handle_wind_h, "dataTime", &data_time);
        codes_handle_delete(handle_wind_h);
    }
    int year = data_date/10000;
    int rest = data_date - 10000*year;
    int month = rest/100;
    rest -= 100*month;
    int day = rest;
    int hour = data_time;
    find_time_coord(year, month, day, hour, 0, 0, 0, t_init);
    if (err != 0)
        ECCERR(err);
    for (int i = 0; i < NUMBER_OF_LEVELS; ++i)
    {
        handle_wind_v = codes_handle_new_from_file(NULL, IN_FILE, PRODUCT_GRIB, &err);
        if (err != 0)
            ECCERR(err);
        if ((retval = codes_get_double_array(handle_wind_v, "values", wind_v, &no_scalars_h)))
            ECCERR(retval);
        codes_handle_delete(handle_wind_v);
        for (int j = 0; j < NUMBER_OF_VECTORS_V; ++j)
            init_state -> wind[j + i*NUMBER_OF_VECTORS_PER_LAYER] = wind_v[j];
    }
    fclose(IN_FILE);
    free(water_vapour_density);
    free(liquid_water_density);
    free(solid_water_density);
    free(liquid_water_temp);
    free(solid_water_temp);
    free(pot_temp);
    free(rho);
    free(wind_h);
    free(wind_v);
    return 0;
}
