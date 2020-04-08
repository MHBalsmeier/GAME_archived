#include "../../../core/src/enum_and_typedefs.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "io.h"
#include "../diagnostics/diagnostics.h"
#include "../r_operators/r_operators.h"
#include "/usr/src/eccodes/include/eccodes.h"
#include "/lib/geos/include/geos.h"
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}

int write_out(State *state_write_out, double t_init, double t_write, char output_foldername[], Grid *grid)
{
    int str_len;
    find_string_length_from_int((int) t_write - t_init, &str_len);
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
    char *SAMPLE_FILE_SCALAR = "../test_generator/grib_files/scalar_field_blueprint_res_id_2.grb2";
    char *SAMPLE_FILE_VECTOR = "../test_generator/grib_files/vector_field_blueprint_res_id_2.grb2";
    FILE *SAMPLE_SCALAR;
    FILE *SAMPLE_VECTOR;
    int err = 0;
    int retval;
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    FILE *OUT_GRIB;
    OUT_GRIB = fopen(OUTPUT_FILE, "w");
    codes_handle *handle_pot_temperature_h = NULL;
    codes_handle *handle_density_h = NULL;
    codes_handle *handle_wind_u_h = NULL;
    codes_handle *handle_wind_v_h = NULL;
    codes_handle *handle_wind_w_h = NULL;
    codes_handle *handle_mslp = NULL;
    double z_height;
    handle_pot_temperature_h = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    handle_density_h = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    handle_mslp = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    SAMPLE_VECTOR = fopen(SAMPLE_FILE_VECTOR, "r");
    handle_wind_u_h = codes_handle_new_from_file(NULL, SAMPLE_VECTOR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_VECTOR);
    SAMPLE_VECTOR = fopen(SAMPLE_FILE_VECTOR, "r");
    handle_wind_v_h = codes_handle_new_from_file(NULL, SAMPLE_VECTOR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_VECTOR);
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    handle_wind_w_h = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    int init_year, init_month, init_day, init_hour, init_minute, init_second, init_microsecond;
    find_hour_from_time_coord(t_init, &init_year, &init_month, &init_day, &init_hour, &init_minute, &init_second, &init_microsecond);
    long data_date = 10000*init_year + 100*init_month + init_day;
    long data_time = init_hour;
    double *pot_temperature_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *rho_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *wind_u_h = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *wind_v_h = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *wind_w_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *mslp = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double pressure_value, mslp_factor, temp_value, temp_mslp, gravity_value, exner_pressure, wind_0, wind_1, wind_u, wind_v;
    double standard_vert_lapse_rate = 0.0065;
    for (int i = 0; i < NUMBER_OF_LAYERS; ++i)
    {
        z_height = grid -> z_scalar[i*NUMBER_OF_SCALARS_H];
        for (int j = 0; j < NUMBER_OF_SCALARS_H; ++j)
        {
            pot_temperature_h[j] = state_write_out -> density_pot_temp[j + i*NUMBER_OF_SCALARS_H]/state_write_out -> density[j + i*NUMBER_OF_SCALARS_H];
            rho_h[j] = state_write_out -> density[j + i*NUMBER_OF_SCALARS_H];
        }
        if (i == NUMBER_OF_LAYERS - 1)
        {
            for (int j = 0; j < NUMBER_OF_SCALARS_H; ++j)
            {
                exner_pressure = pow(R_D*state_write_out -> density_pot_temp[j + i*NUMBER_OF_SCALARS_H]/P_0, R_D/C_V);
                temp_value = exner_pressure*state_write_out -> density_pot_temp[j + i*NUMBER_OF_SCALARS_H]/state_write_out -> density[j + i*NUMBER_OF_SCALARS_H];
                pressure_value = state_write_out -> density[j + i*NUMBER_OF_SCALARS_H]*R_D*temp_value;
                temp_mslp = temp_value + standard_vert_lapse_rate*grid -> z_scalar[j + i*NUMBER_OF_SCALARS_H];
                gravity_value = 9.80616;
                mslp_factor = pow(1 - (temp_mslp - temp_value)/temp_mslp, gravity_value/(R_D*standard_vert_lapse_rate));
                mslp[j] = pressure_value/mslp_factor;
            }
            if ((retval = codes_set_long(handle_mslp, "discipline", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "centre", 255)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "significanceOfReferenceTime", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "productionStatusOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "typeOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "indicatorOfUnitOfTimeRange", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "stepUnits", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "dataDate", data_date)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "dataTime", data_time)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "forecastTime", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "stepRange", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "typeOfGeneratingProcess", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "parameterCategory", 3)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "parameterNumber", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "typeOfFirstFixedSurface", 101)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "scaleFactorOfFirstFixedSurface", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "scaledValueOfFirstFixedSurface", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "level", i + 1)))
                ECCERR(retval);
            if ((retval = codes_set_double_array(handle_mslp, "values", mslp, NUMBER_OF_SCALARS_H)))
                ECCERR(retval);
            if ((retval = codes_write_message(handle_mslp, OUTPUT_FILE, "a")))
                ECCERR(retval);
        }
        if ((retval = codes_set_long(handle_pot_temperature_h, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "stepUnits", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "dataDate", data_date)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "dataTime", data_time)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "forecastTime", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "stepRange", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "parameterCategory", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "parameterNumber", 2)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "scaledValueOfFirstFixedSurface", z_height)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_pot_temperature_h, "values", pot_temperature_h, NUMBER_OF_SCALARS_H)))
            ECCERR(retval);
        if (i == 0)
        {
            if ((retval = codes_write_message(handle_pot_temperature_h, OUTPUT_FILE, "w")))
                ECCERR(retval);
        }
        else
        {
            if ((retval = codes_write_message(handle_pot_temperature_h, OUTPUT_FILE, "a")))
                ECCERR(retval);
        }
        if ((retval = codes_set_long(handle_density_h, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "stepUnits", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "dataDate", data_date)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "dataTime", data_time)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "forecastTime", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "stepRange", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "parameterCategory", 3)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "parameterNumber", 10)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "scaledValueOfFirstFixedSurface", z_height)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_density_h, "values", rho_h, NUMBER_OF_SCALARS_H)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_density_h, OUTPUT_FILE, "a")))
            ECCERR(retval);
        for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
        {
            wind_0 = state_write_out -> wind[j + i*NUMBER_OF_VECTORS_H + (i + 1)*NUMBER_OF_VECTORS_V];
            retval = recov_hor_par_pri(state_write_out -> wind, i, j, &wind_1, grid);
            retval = passive_turn(wind_0, wind_1, -grid -> direction[j], &wind_u, &wind_v);
            wind_u_h[j] = wind_u;
            wind_v_h[j] = wind_v;
        }
        if ((retval = codes_set_long(handle_wind_u_h, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "stepUnits", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "dataDate", data_date)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "dataTime", data_time)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "forecastTime", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "stepRange", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "parameterCategory", 2)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "parameterNumber", 2)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "scaledValueOfFirstFixedSurface", z_height)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_wind_u_h, "values", wind_u_h, NUMBER_OF_VECTORS_H)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_wind_u_h, OUTPUT_FILE, "a")))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "stepUnits", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "dataDate", data_date)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "dataTime", data_time)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "forecastTime", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "stepRange", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "parameterCategory", 2)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "parameterNumber", 3)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "scaledValueOfFirstFixedSurface", z_height)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_wind_v_h, "values", wind_v_h, NUMBER_OF_VECTORS_H)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_wind_v_h, OUTPUT_FILE, "a")))
            ECCERR(retval);
    }
    free(mslp);
    codes_handle_delete(handle_mslp);
    codes_handle_delete(handle_pot_temperature_h);
    codes_handle_delete(handle_density_h);
    codes_handle_delete(handle_wind_u_h);
    codes_handle_delete(handle_wind_v_h);
    SAMPLE_SCALAR = fopen(SAMPLE_FILE_SCALAR, "r");
    handle_wind_v_h = codes_handle_new_from_file(NULL, SAMPLE_SCALAR, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE_SCALAR);
    for (int i = 0; i < NUMBER_OF_LEVELS; ++i)
    {
        z_height = grid -> z_vector[i*NUMBER_OF_VECTORS_PER_LAYER];
        for (int j = 0; j < NUMBER_OF_SCALARS_H; j++)
            wind_w_h[j] = state_write_out -> wind[j + i*(NUMBER_OF_SCALARS_H + NUMBER_OF_VECTORS_H)];
        if ((retval = codes_set_long(handle_wind_w_h, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "stepUnits", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "dataDate", data_date)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "dataTime", data_time)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "forecastTime", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "stepRange", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "parameterCategory", 2)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "parameterNumber", 9)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "scaledValueOfFirstFixedSurface", z_height)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_wind_w_h, "values", wind_v_h, NUMBER_OF_SCALARS_H)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_wind_w_h, OUTPUT_FILE, "a")))
            ECCERR(retval);
    }
    codes_handle_delete(handle_wind_w_h);
    free(pot_temperature_h);
    free(rho_h);
    free(wind_u_h);
    free(wind_v_h);
    free(wind_w_h);
    fclose(OUT_GRIB);
    return 0;
}

int find_string_length_from_int(int input, int *answer)
{
    *answer = 1;
    for (int i = 1; i < 10; ++i)
    {
        if (input >= pow(10, i))
            ++*answer;
    }
    return 0;
}
