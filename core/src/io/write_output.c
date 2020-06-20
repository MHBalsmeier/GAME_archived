#include "../../../core/src/enum_and_typedefs.h"
#include <stdio.h>
#include <string.h>
#include "io.h"
#include "../diagnostics/diagnostics.h"
#include "../spatial_operators/spatial_operators.h"
#include "eccodes.h"
#include "geos95.h"
#include "atmostracers.h"
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
    char *SAMPLE_FILENAME = "./input/grib_template.grb2";
    FILE *SAMPLE_FILE;
    int err = 0;
    int retval;
    FILE *OUT_GRIB;
    OUT_GRIB = fopen(OUTPUT_FILE, "w+");
    int init_year, init_month, init_day, init_hour, init_minute, init_second, init_microsecond;
    find_hour_from_time_coord(t_init, &init_year, &init_month, &init_day, &init_hour, &init_minute, &init_second, &init_microsecond);
    long data_date = 10000*init_year + 100*init_month + init_day;
    long data_time = init_hour;
    double *pot_temperature_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *rho_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *wind_u_h = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *wind_v_h = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *wind_w_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *wind_u_10m = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *wind_v_10m = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *mslp = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *surface_p = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *t2 = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *tcdc = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *rprate = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *sprate = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double pressure_value, mslp_factor, surface_p_factor, temp_lowest_layer, temp_mslp, temp_surface, exner_pressure, wind_0, wind_1, wind_u, wind_v, delta_z_temp, temp_gradient, temp_upper, temp_lower;
    double standard_vert_lapse_rate = 0.0065;
    double pot_temp, condensates_density_sum, density_d_micro_value, density_v_micro_value;
    long unsigned length = 4;
    codes_handle *handle_pot_temperature_h = NULL;
    codes_handle *handle_density_h = NULL;
    codes_handle *handle_wind_u_h = NULL;
    codes_handle *handle_wind_v_h = NULL;
    codes_handle *handle_wind_w_h = NULL;
    codes_handle *handle_wind_u_10m = NULL;
    codes_handle *handle_wind_v_10m = NULL;
    codes_handle *handle_mslp = NULL;
    codes_handle *handle_surface_p = NULL;
    codes_handle *handle_t2 = NULL;
    codes_handle *handle_tcdc = NULL;
    codes_handle *handle_rprate = NULL;
    codes_handle *handle_sprate = NULL;
    for (int i = 0; i < NUMBER_OF_LAYERS; ++i)
    {
        for (int j = 0; j < NUMBER_OF_SCALARS_H; ++j)
        {
        	rho_h[j] = state_write_out -> density_dry[j + i*NUMBER_OF_SCALARS_H];
        	condensates_density_sum = calc_condensates_density_sum(i, j, state_write_out -> tracer_densities);
			pot_temperature_h[j] = pot_temp_diagnostics_single_value(state_write_out -> entropy_gas[j + i*NUMBER_OF_SCALARS_H], rho_h[j], state_write_out -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + j + i*NUMBER_OF_SCALARS_H], condensates_density_sum);
        }
        if (i == NUMBER_OF_LAYERS - 1)
        {
            for (int j = 0; j < NUMBER_OF_SCALARS_H; ++j)
            {
            	condensates_density_sum = calc_condensates_density_sum(i, j, state_write_out -> tracer_densities);
            	pot_temp = pot_temp_diagnostics_single_value(state_write_out -> entropy_gas[j + i*NUMBER_OF_SCALARS_H], state_write_out -> density_dry[j + i*NUMBER_OF_SCALARS_H], state_write_out -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + j + i*NUMBER_OF_SCALARS_H], condensates_density_sum);
            	density_d_micro_value = calc_micro_density(state_write_out -> density_dry[j + i*NUMBER_OF_SCALARS_H], condensates_density_sum);
            	density_v_micro_value = calc_micro_density(state_write_out -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + j + i*NUMBER_OF_SCALARS_H], condensates_density_sum);
                exner_pressure = exner_pressure_diagnostics_single_value(density_d_micro_value, density_v_micro_value, pot_temp);
                temp_lowest_layer = temperature_diagnostics_single_value(exner_pressure, pot_temp);
                pressure_value = state_write_out -> density_dry[j + i*NUMBER_OF_SCALARS_H]*R_D*temp_lowest_layer;
                temp_mslp = temp_lowest_layer + standard_vert_lapse_rate*grid -> z_scalar[j + i*NUMBER_OF_SCALARS_H];
                mslp_factor = pow(1 - (temp_mslp - temp_lowest_layer)/temp_mslp, -grid -> gravity[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_PER_LAYER + j]/(R_D*standard_vert_lapse_rate));
                mslp[j] = pressure_value/mslp_factor;
				temp_surface = temp_lowest_layer + standard_vert_lapse_rate*(grid -> z_scalar[j + i*NUMBER_OF_SCALARS_H] - grid -> z_vector[NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_V + j]);
                surface_p_factor = pow(1 - (temp_surface - temp_lowest_layer)/temp_surface, -grid -> gravity[NUMBER_OF_LAYERS*NUMBER_OF_VECTORS_PER_LAYER + j]/(R_D*standard_vert_lapse_rate));
				surface_p[j] = pressure_value/surface_p_factor;
                delta_z_temp = 2 - grid -> z_scalar[j + i*NUMBER_OF_SCALARS_H];
            	condensates_density_sum = calc_condensates_density_sum(i - 1, j, state_write_out -> tracer_densities);
                pot_temp = pot_temp_diagnostics_single_value(state_write_out -> entropy_gas[j + (i - 1)*NUMBER_OF_SCALARS_H], state_write_out -> density_dry[j + (i - 1)*NUMBER_OF_SCALARS_H], state_write_out -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + j + (i - 1)*NUMBER_OF_SCALARS_H], condensates_density_sum);
            	density_d_micro_value = calc_micro_density(state_write_out -> density_dry[j + (i - 1)*NUMBER_OF_SCALARS_H], condensates_density_sum);
            	density_v_micro_value = calc_micro_density(state_write_out -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + j + (i - 1)*NUMBER_OF_SCALARS_H], condensates_density_sum);
                exner_pressure = exner_pressure_diagnostics_single_value(density_d_micro_value, density_v_micro_value, pot_temp);
                temp_upper = temperature_diagnostics_single_value(exner_pressure, pot_temp);
            	condensates_density_sum = calc_condensates_density_sum(i, j, state_write_out -> tracer_densities);
                pot_temp = pot_temp_diagnostics_single_value(state_write_out -> entropy_gas[j + i*NUMBER_OF_SCALARS_H], state_write_out -> density_dry[j + i*NUMBER_OF_SCALARS_H], state_write_out -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i*NUMBER_OF_SCALARS_H + j], condensates_density_sum);
            	density_d_micro_value = calc_micro_density(state_write_out -> density_dry[j + i*NUMBER_OF_SCALARS_H], condensates_density_sum);
            	density_v_micro_value = calc_micro_density(state_write_out -> tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + j + i*NUMBER_OF_SCALARS_H], condensates_density_sum);
                exner_pressure = exner_pressure_diagnostics_single_value(density_d_micro_value, density_v_micro_value, pot_temp);
                temp_lower = temperature_diagnostics_single_value(exner_pressure, pot_temp);
                temp_gradient = (temp_upper - temp_lower)/(grid -> z_scalar[j + (i - 1)*NUMBER_OF_SCALARS_H] - grid -> z_scalar[j + i*NUMBER_OF_SCALARS_H]);
                t2[j] = temp_lowest_layer + delta_z_temp*temp_gradient;
                sprate[j] = 0;
                rprate[j] = 0;
                tcdc[j] = 0;
                for (int k = 0; k < NUMBER_OF_CONDENSATED_TRACERS; ++k)
                {
                    for (int l = 0; l < NUMBER_OF_LAYERS; ++l)
                    {
                        if (state_write_out -> tracer_densities[k*NUMBER_OF_SCALARS + l*NUMBER_OF_SCALARS_H + j] > 0)
                            tcdc[j] = 100*1;
                    }
                }
                for (int k = 0; k < NUMBER_OF_SOLID_TRACERS; ++k)
                    sprate[j] += fmax(ret_sink_velocity(k, 0, 0.001)*state_write_out -> tracer_densities[k*NUMBER_OF_SCALARS + i*NUMBER_OF_SCALARS_H + j], 0);
                for (int k = NUMBER_OF_SOLID_TRACERS; k < NUMBER_OF_CONDENSATED_TRACERS; ++k)
                    rprate[j] += fmax(ret_sink_velocity(k, 0, 0.001)*state_write_out -> tracer_densities[k*NUMBER_OF_SCALARS + i*NUMBER_OF_SCALARS_H + j], 0);
            }
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_mslp = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
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
            if ((retval = codes_set_long(handle_mslp, "typeOfFirstFixedSurface", 102)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "scaledValueOfFirstFixedSurface", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "scaleFactorOfFirstFixedSurface", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_mslp, "level", 0)))
                ECCERR(retval);
            if ((retval = codes_set_double_array(handle_mslp, "values", mslp, NUMBER_OF_SCALARS_H)))
                ECCERR(retval);
			codes_handle_delete(handle_mslp);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_surface_p = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
            if ((retval = codes_write_message(handle_surface_p, OUTPUT_FILE, "a")))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "discipline", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "centre", 255)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "significanceOfReferenceTime", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "productionStatusOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "typeOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "indicatorOfUnitOfTimeRange", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "stepUnits", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "dataDate", data_date)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "dataTime", data_time)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "forecastTime", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "stepRange", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "typeOfGeneratingProcess", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "parameterCategory", 3)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "parameterNumber", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "typeOfFirstFixedSurface", 103)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "scaledValueOfFirstFixedSurface", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "scaleFactorOfFirstFixedSurface", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "level", 0)))
                ECCERR(retval);
            if ((retval = codes_set_double_array(handle_surface_p, "values", surface_p, NUMBER_OF_SCALARS_H)))
                ECCERR(retval);
            if ((retval = codes_write_message(handle_surface_p, OUTPUT_FILE, "a")))
                ECCERR(retval);
			codes_handle_delete(handle_surface_p);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_t2 = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
            if ((retval = codes_set_long(handle_t2, "discipline", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "centre", 255)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "significanceOfReferenceTime", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "productionStatusOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "typeOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "indicatorOfUnitOfTimeRange", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "stepUnits", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "dataDate", data_date)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "dataTime", data_time)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "forecastTime", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "stepRange", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "typeOfGeneratingProcess", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "parameterCategory", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "parameterNumber", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "typeOfFirstFixedSurface", 103)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "scaledValueOfFirstFixedSurface", 2)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "scaleFactorOfFirstFixedSurface", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_t2, "level", 2)))
                ECCERR(retval);
            if ((retval = codes_set_double_array(handle_t2, "values", t2, NUMBER_OF_SCALARS_H)))
                ECCERR(retval);
            if ((retval = codes_write_message(handle_t2, OUTPUT_FILE, "a")))
                ECCERR(retval);
			codes_handle_delete(handle_t2);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_rprate = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
            if ((retval = codes_set_long(handle_rprate, "discipline", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "centre", 255)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "significanceOfReferenceTime", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "productionStatusOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "typeOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "indicatorOfUnitOfTimeRange", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "stepUnits", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "dataDate", data_date)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "dataTime", data_time)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "forecastTime", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "stepRange", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "typeOfGeneratingProcess", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "parameterCategory", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "parameterNumber", 65)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "typeOfFirstFixedSurface", 103)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "scaledValueOfFirstFixedSurface", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "scaleFactorOfFirstFixedSurface", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_rprate, "level", 0)))
                ECCERR(retval);
            if ((retval = codes_set_double_array(handle_rprate, "values", rprate, NUMBER_OF_SCALARS_H)))
                ECCERR(retval);
            if ((retval = codes_write_message(handle_rprate, OUTPUT_FILE, "a")))
                ECCERR(retval);
			codes_handle_delete(handle_rprate);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_sprate = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
            if ((retval = codes_set_long(handle_sprate, "discipline", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "centre", 255)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "significanceOfReferenceTime", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "productionStatusOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "typeOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "indicatorOfUnitOfTimeRange", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "stepUnits", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "dataDate", data_date)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "dataTime", data_time)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "forecastTime", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "stepRange", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "typeOfGeneratingProcess", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "parameterCategory", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "parameterNumber", 66)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "typeOfFirstFixedSurface", 103)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "scaledValueOfFirstFixedSurface", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "scaleFactorOfFirstFixedSurface", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_sprate, "level", 0)))
                ECCERR(retval);
            if ((retval = codes_set_double_array(handle_sprate, "values", sprate, NUMBER_OF_SCALARS_H)))
                ECCERR(retval);
            if ((retval = codes_write_message(handle_sprate, OUTPUT_FILE, "a")))
                ECCERR(retval);
			codes_handle_delete(handle_sprate);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_tcdc = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
            if ((retval = codes_set_long(handle_tcdc, "discipline", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "centre", 255)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "significanceOfReferenceTime", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "productionStatusOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "typeOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "indicatorOfUnitOfTimeRange", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "stepUnits", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "dataDate", data_date)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "dataTime", data_time)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "forecastTime", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "stepRange", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "typeOfGeneratingProcess", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "parameterCategory", 6)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "parameterNumber", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "typeOfFirstFixedSurface", 103)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "scaledValueOfFirstFixedSurface", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "scaleFactorOfFirstFixedSurface", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_tcdc, "level", 0)))
                ECCERR(retval);
            if ((retval = codes_set_string(handle_tcdc, "shortName", "tcc", &length)))
                ECCERR(retval);
            if ((retval = codes_set_double_array(handle_tcdc, "values", tcdc, NUMBER_OF_SCALARS_H)))
                ECCERR(retval);
            if ((retval = codes_write_message(handle_tcdc, OUTPUT_FILE, "a")))
                ECCERR(retval);
			codes_handle_delete(handle_tcdc);
        }
    	SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_pot_temperature_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
		    ECCERR(err);
		fclose(SAMPLE_FILE);
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
        if ((retval = codes_set_long(handle_pot_temperature_h, "typeOfFirstFixedSurface", 26)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "scaledValueOfFirstFixedSurface", i)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature_h, "scaleFactorOfFirstFixedSurface", 1)))
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
		codes_handle_delete(handle_pot_temperature_h);
		SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_density_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
		    ECCERR(err);
		fclose(SAMPLE_FILE);
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
        if ((retval = codes_set_long(handle_density_h, "typeOfFirstFixedSurface", 26)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "scaledValueOfFirstFixedSurface", i)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density_h, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_density_h, "values", rho_h, NUMBER_OF_SCALARS_H)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_density_h, OUTPUT_FILE, "a")))
            ECCERR(retval);
		codes_handle_delete(handle_density_h);
        for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
        {
            wind_0 = state_write_out -> velocity_gas[j + i*NUMBER_OF_VECTORS_H + (i + 1)*NUMBER_OF_VECTORS_V];
            retval = recov_hor_par_pri(state_write_out -> velocity_gas, i, j, &wind_1, grid);
            retval = passive_turn(wind_0, wind_1, -grid -> direction[j], &wind_u, &wind_v);
            wind_u_h[j] = wind_u;
            wind_v_h[j] = wind_v;
        }
		SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_wind_u_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
		    ECCERR(err);
		fclose(SAMPLE_FILE);
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
        if ((retval = codes_set_long(handle_wind_u_h, "typeOfFirstFixedSurface", 26)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "scaledValueOfFirstFixedSurface", i)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_u_h, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_wind_u_h, "values", wind_u_h, NUMBER_OF_VECTORS_H)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_wind_u_h, OUTPUT_FILE, "a")))
            ECCERR(retval);
		codes_handle_delete(handle_wind_u_h);
		SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_wind_v_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
		    ECCERR(err);
		fclose(SAMPLE_FILE);
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
        if ((retval = codes_set_long(handle_wind_v_h, "typeOfFirstFixedSurface", 26)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "scaledValueOfFirstFixedSurface", i)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v_h, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_wind_v_h, "values", wind_v_h, NUMBER_OF_VECTORS_H)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_wind_v_h, OUTPUT_FILE, "a")))
            ECCERR(retval);
		codes_handle_delete(handle_wind_v_h);
        if (i == NUMBER_OF_LAYERS - 1)
        {
            for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
            {
                wind_u_10m[j] = 0.667*wind_u_h[j];
                wind_v_10m[j] = 0.667*wind_v_h[j];
            }
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_wind_u_10m = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
            if ((retval = codes_set_long(handle_wind_u_10m, "discipline", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "centre", 255)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "significanceOfReferenceTime", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "productionStatusOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "typeOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "indicatorOfUnitOfTimeRange", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "stepUnits", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "dataDate", data_date)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "dataTime", data_time)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "forecastTime", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "stepRange", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "typeOfGeneratingProcess", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "parameterCategory", 2)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "parameterNumber", 2)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "typeOfFirstFixedSurface", 103)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "scaledValueOfFirstFixedSurface", 10)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "scaleFactorOfFirstFixedSurface", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_u_10m, "level", 10)))
                ECCERR(retval);
            if ((retval = codes_set_double_array(handle_wind_u_10m, "values", wind_u_10m, NUMBER_OF_VECTORS_H)))
                ECCERR(retval);
            if ((retval = codes_write_message(handle_wind_u_10m, OUTPUT_FILE, "a")))
                ECCERR(retval);
			codes_handle_delete(handle_wind_u_10m);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_wind_v_10m = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
            if ((retval = codes_set_long(handle_wind_v_10m, "discipline", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "centre", 255)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "significanceOfReferenceTime", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "productionStatusOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "typeOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "indicatorOfUnitOfTimeRange", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "stepUnits", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "dataDate", data_date)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "dataTime", data_time)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "forecastTime", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "stepRange", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "typeOfGeneratingProcess", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "parameterCategory", 2)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "parameterNumber", 3)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "typeOfFirstFixedSurface", 103)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "scaledValueOfFirstFixedSurface", 10)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "scaleFactorOfFirstFixedSurface", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_wind_v_10m, "level", 10)))
                ECCERR(retval);
            if ((retval = codes_set_double_array(handle_wind_v_10m, "values", wind_v_10m, NUMBER_OF_VECTORS_H)))
                ECCERR(retval);
            if ((retval = codes_write_message(handle_wind_v_10m, OUTPUT_FILE, "a")))
                ECCERR(retval);
			codes_handle_delete(handle_wind_v_10m);
        }
    }
    free(t2);
    free(mslp);
    free(surface_p);
    free(rprate);
    free(sprate);
    free(tcdc);
    free(pot_temperature_h);
    free(rho_h);
    free(wind_u_h);
    free(wind_v_h);
    free(wind_u_10m);
    free(wind_v_10m);
	SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
	handle_wind_w_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
	    ECCERR(err);
	fclose(SAMPLE_FILE);
    for (int i = 0; i < NUMBER_OF_LEVELS; ++i)
    {
        for (int j = 0; j < NUMBER_OF_VECTORS_V; j++)
            wind_w_h[j] = state_write_out -> velocity_gas[j + i*NUMBER_OF_VECTORS_PER_LAYER];
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
        if ((retval = codes_set_long(handle_wind_w_h, "typeOfFirstFixedSurface", 26)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "scaledValueOfFirstFixedSurface", i)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_w_h, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_wind_w_h, "values", wind_v_h, NUMBER_OF_VECTORS_V)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_wind_w_h, OUTPUT_FILE, "a")))
            ECCERR(retval);
    }
    codes_handle_delete(handle_wind_w_h);
    free(wind_w_h);
    fclose(OUT_GRIB);
    return 0;
}

int write_out_integral(State *state_write_out, double t_write, char output_foldername[], Grid *grid, Dualgrid *dualgrid, int integral_id)
{
	/*
	integral_id:
	0: dry mass
	1: entropy
	2: energy
	*/
    double global_integral = 0;
    FILE *global_integral_file;
    int INTEGRAL_FILE_LENGTH = 200;
    char *INTEGRAL_FILE_PRE = malloc((INTEGRAL_FILE_LENGTH + 1)*sizeof(char));
    if (integral_id == 0)
   		sprintf(INTEGRAL_FILE_PRE, "%s/%s", output_foldername, "dry_mass");
    if (integral_id == 1)
   		sprintf(INTEGRAL_FILE_PRE, "%s/%s", output_foldername, "entropy");
    if (integral_id == 2)
   		sprintf(INTEGRAL_FILE_PRE, "%s/%s", output_foldername, "energy");
    INTEGRAL_FILE_LENGTH = strlen(INTEGRAL_FILE_PRE);
    char *INTEGRAL_FILE = malloc((INTEGRAL_FILE_LENGTH + 1)*sizeof(char));
    sprintf(INTEGRAL_FILE, "%s", INTEGRAL_FILE_PRE);
    free(INTEGRAL_FILE_PRE);
    int retval;
    if (integral_id == 0)
    {
    	global_integral_file = fopen(INTEGRAL_FILE, "a");
    	global_scalar_integrator(state_write_out -> density_dry, grid, &global_integral);
    	fprintf(global_integral_file, "%lf\t%lf\n", t_write, global_integral);
    	fclose(global_integral_file);
    }
    if (integral_id == 1)
    {
    	global_integral_file = fopen(INTEGRAL_FILE, "a");
    	global_scalar_integrator(state_write_out -> entropy_gas, grid, &global_integral);
    	fprintf(global_integral_file, "%lf\t%lf\n", t_write, global_integral);
    	fclose(global_integral_file);
    }
    if (integral_id == 2)
    {
    	double kinetic_integral, potential_integral, internal_integral;
    	global_integral_file = fopen(INTEGRAL_FILE, "a");
    	Scalar_field *e_kin_density = malloc(sizeof(Scalar_field));
    	retval = kinetic_energy(state_write_out -> velocity_gas, *e_kin_density, grid);
    	if (retval != 0)
    	{
    		printf("Error in kinetic_energy called from write_output, position 0. Answer is %d.\n", retval);
    		exit(1);
    	}
    	retval = scalar_times_scalar(state_write_out -> density_dry, *e_kin_density, *e_kin_density);
    	global_scalar_integrator(*e_kin_density, grid, &kinetic_integral);
    	free(e_kin_density);
    	Scalar_field *pot_energy_density = malloc(sizeof(Scalar_field));
    	retval = scalar_times_scalar(state_write_out -> density_dry, grid -> gravity_potential, *pot_energy_density);
    	global_scalar_integrator(*pot_energy_density, grid, &potential_integral);
    	free(pot_energy_density);
    	Scalar_field *int_energy_density = malloc(sizeof(Scalar_field));
    	retval = temperature_diagnostics(state_write_out -> entropy_gas, state_write_out -> density_dry, state_write_out -> tracer_densities, *int_energy_density);
    	retval = scalar_times_scalar(state_write_out -> density_dry, *int_energy_density, *int_energy_density);
    	global_scalar_integrator(*int_energy_density, grid, &internal_integral);
    	fprintf(global_integral_file, "%lf\t%lf\t%lf\t%lf\n", t_write, kinetic_integral, potential_integral, C_D_V*internal_integral);
    	free(int_energy_density);
    	fclose(global_integral_file);
    }
    free(INTEGRAL_FILE);
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

