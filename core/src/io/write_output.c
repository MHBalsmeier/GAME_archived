#include "/home/max/my_code/game/core/src/enum_and_typedefs.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "io.h"
#include "../diagnostics/diagnostics.h"
#include "/home/max/custom_builds/eccodes/include/eccodes.h"
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}

int write_out(State *state_write_out, double t_init, double t_write, char output_foldername[])
{
    const double ATMOS_HEIGHT = SCALE_HEIGHT*log(1 + NUMBER_OF_LAYERS);
    int str_len;
    find_string_length_from_int((int) t_write, &str_len);
    char add_time_str[str_len];
    char pre_string[] = "init+";
    char append_string[] = "s.grb2";
    sprintf(add_time_str, "%lf", t_write);
    int out_file_length = strlen(output_foldername) + 1 + strlen(pre_string) + str_len + strlen(append_string);
    char OUTPUT_FILE[out_file_length];
    for (int i = 0; i < out_file_length + 1; i++)
    {
        if (i < strlen(output_foldername))
            OUTPUT_FILE[i] = output_foldername[i];
        if (i == strlen(output_foldername))
            OUTPUT_FILE[i] = '/';
        if (i >= strlen(output_foldername) + 1 && i < strlen(output_foldername) + 1 + strlen(pre_string))
            OUTPUT_FILE[i] = pre_string[i - (strlen(output_foldername) + 1)];
        if (i >= strlen(output_foldername) + 1 + strlen(pre_string) && i < out_file_length - strlen(append_string))
            OUTPUT_FILE[i] = add_time_str[i - (strlen(output_foldername) + 1 + strlen(pre_string))];
        if (i >= out_file_length - strlen(append_string))
            OUTPUT_FILE[i] = append_string[i - (out_file_length - strlen(append_string))];
    }
    int dimids_scalar[2];
    int scalar_index;
    const int START_SECTION = 4;
    char *SAMPLE_FILE_SCALAR = "/home/max/my_code/game/test_generator/grib_files/scalar_field_blueprint_res_id_2.grb2";
    char *SAMPLE_FILE_VECTOR = "/home/max/my_code/game/test_generator/grib_files/vector_field_blueprint_res_id_2.grb2";
    FILE *SAMPLE_SCALAR;
    FILE *SAMPLE_VECTOR;
    int err = 0;
    int retval;
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    FILE *OUT_GRIB;
    OUT_GRIB = fopen(OUTPUT_FILE, "w");
    codes_handle *handle_pot_temperature_h = NULL;
    codes_handle *handle_density_h = NULL;
    codes_handle *handle_wind_h_h = NULL;
    codes_handle *handle_wind_v_h = NULL;
    double *pot_temperature_h, *rho_h, *wind_h_h, *wind_v_h;
    pot_temperature_h = malloc(sizeof(double)*NUMBER_OF_SCALARS_H);
    rho_h = malloc(sizeof(double)*NUMBER_OF_SCALARS_H);
    wind_h_h = malloc(sizeof(double)*NUMBER_OF_VECTORS_H);
    wind_v_h = malloc(sizeof(double)*NUMBER_OF_SCALARS_H);
    double sigma, z_height;
    handle_pot_temperature_h = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    handle_density_h = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    SAMPLE_VECTOR = fopen(SAMPLE_FILE_VECTOR, "r");
    handle_wind_h_h = codes_handle_new_from_file(NULL, SAMPLE_VECTOR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_VECTOR);
    int init_year, init_month, init_day, init_hour, init_minute, init_second, init_microsecond;
    find_hour_from_time_coord(t_init, &init_year, &init_month, &init_day, &init_hour, &init_minute, &init_second, &init_microsecond);
    long data_date = 10000*init_year + 100*init_month + init_day;
    long data_time = init_hour;
    for (int i = 0; i < NUMBER_OF_LAYERS; ++i)
    {
        for (int j = 0; j < NUMBER_OF_SCALARS_H; ++j)
        {
            pot_temperature_h[j] = state_write_out -> density_pot_temp[j + i*NUMBER_OF_SCALARS_H]/state_write_out -> density[j + i*NUMBER_OF_SCALARS_H];
            rho_h[j] = state_write_out -> density[j + i*NUMBER_OF_SCALARS_H];
        }
        if (retval = codes_set_long(handle_pot_temperature_h, "discipline", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "centre", 255))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "significanceOfReferenceTime", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "productionStatusOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "typeOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "indicatorOfUnitOfTimeRange", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "stepUnits", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "dataDate", data_date))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "dataTime", data_time))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "forecastTime", t_write - t_init))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "stepRange", t_write - t_init))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "typeOfGeneratingProcess", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "parameterCategory", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "parameterNumber", 2))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "typeOfFirstFixedSurface", 102))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "scaleFactorOfFirstFixedSurface", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "scaledValueOfFirstFixedSurface", z_height))
            ECCERR(retval);
        if (retval = codes_set_long(handle_pot_temperature_h, "level", i))
            ECCERR(retval);
        if (retval = codes_set_double_array(handle_pot_temperature_h, "values", pot_temperature_h, NUMBER_OF_SCALARS_H))
            ECCERR(retval);
        if (i == 0)
            codes_write_message(handle_pot_temperature_h, OUTPUT_FILE, "w");
        else
            codes_write_message(handle_pot_temperature_h, OUTPUT_FILE, "a");
        if (retval = codes_set_long(handle_density_h, "discipline", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "centre", 255))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "significanceOfReferenceTime", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "productionStatusOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "typeOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "indicatorOfUnitOfTimeRange", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "stepUnits", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "dataDate", data_date))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "dataTime", data_time))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "forecastTime", t_write - t_init))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "stepRange", t_write - t_init))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "typeOfGeneratingProcess", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "parameterCategory", 3))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "parameterNumber", 10))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "typeOfFirstFixedSurface", 102))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "scaleFactorOfFirstFixedSurface", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "scaledValueOfFirstFixedSurface", z_height))
            ECCERR(retval);
        if (retval = codes_set_long(handle_density_h, "level", i))
            ECCERR(retval);
        if (retval = codes_set_double_array(handle_density_h, "values", rho_h, NUMBER_OF_SCALARS_H))
            ECCERR(retval);
        codes_write_message(handle_density_h, OUTPUT_FILE, "a");
        for (int j = 0; j < NUMBER_OF_VECTORS_H; j++)
            wind_h_h[j] = state_write_out -> wind[j + i*NUMBER_OF_VECTORS_H + (i + 1)*NUMBER_OF_SCALARS_H];
        if (retval = codes_set_long(handle_wind_h_h, "discipline", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "centre", 255))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "significanceOfReferenceTime", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "productionStatusOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "typeOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "indicatorOfUnitOfTimeRange", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "stepUnits", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "dataDate", data_date))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "dataTime", data_time))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "forecastTime", t_write - t_init))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "stepRange", t_write - t_init))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "typeOfGeneratingProcess", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "parameterCategory", 2))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "parameterNumber", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "typeOfFirstFixedSurface", 102))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "scaleFactorOfFirstFixedSurface", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "scaledValueOfFirstFixedSurface", z_height))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_h_h, "level", i))
            ECCERR(retval);
        if (retval = codes_set_double_array(handle_wind_h_h, "values", wind_h_h, NUMBER_OF_VECTORS_H))
            ECCERR(retval);
        codes_write_message(handle_wind_h_h, OUTPUT_FILE, "a");
    }
    codes_handle_delete(handle_pot_temperature_h);
    codes_handle_delete(handle_density_h);
    codes_handle_delete(handle_wind_h_h);
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    handle_wind_v_h = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    for (int i = 0; i < NUMBER_OF_LEVELS; i++)
    {
        sigma = (SCALE_HEIGHT/ATMOS_HEIGHT)*log((1.0 + NUMBER_OF_LAYERS)/(i + 1));
        z_height = ATMOS_HEIGHT*sigma;
        for (int j = 0; j < NUMBER_OF_SCALARS_H; j++)
            wind_v_h[j] = state_write_out -> wind[j + i*(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H)];
        if (retval = codes_set_long(handle_wind_v_h, "discipline", 0))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "centre", 255))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "significanceOfReferenceTime", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "productionStatusOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "typeOfProcessedData", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "indicatorOfUnitOfTimeRange", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "stepUnits", 13))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "dataDate", data_date))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "dataTime", data_time))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "forecastTime", t_write - t_init))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "stepRange", t_write - t_init))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "typeOfGeneratingProcess", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "parameterCategory", 2))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "parameterNumber", 9))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "typeOfFirstFixedSurface", 102))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "scaleFactorOfFirstFixedSurface", 1))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "scaledValueOfFirstFixedSurface", z_height))
            ECCERR(retval);
        if (retval = codes_set_long(handle_wind_v_h, "level", i))
            ECCERR(retval);
        if (retval = codes_set_double_array(handle_wind_v_h, "values", wind_v_h, NUMBER_OF_SCALARS_H))
            ECCERR(retval);
        codes_write_message(handle_wind_v_h, OUTPUT_FILE, "a");
    }
    codes_handle_delete(handle_wind_v_h);
    free(pot_temperature_h);
    free(rho_h);
    free(wind_h_h);
    free(wind_v_h);
    fclose(OUT_GRIB);
    return 0;
}

int find_string_length_from_int(int input, int *answer)
{
    *answer = 1;
    for (int i = 1; i < 10; ++i)
    {
        if (input > pow(10, i))
            ++*answer;
    }
    return 0;
}
