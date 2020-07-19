/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
Here, the output is written to grib files and integrals are written to text files if configured that way.
*/

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

double calc_std_dev(double [], int);

int write_out(State *state_write_out, double wind_h_lowest_layer_array[], int min_no_of_output_steps, double t_init, double t_write, char output_directory[], Grid *grid, Dualgrid *dualgrid)
{
    int init_year, init_month, init_day, init_hour, init_minute, init_second, init_microsecond;
    find_hour_from_time_coord(t_init, &init_year, &init_month, &init_day, &init_hour, &init_minute, &init_second, &init_microsecond);
    long data_date = 10000*init_year + 100*init_month + init_day;
    long data_time = init_hour;
    int OUTPUT_FILE_LENGTH = 300;
    char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE_PRE, "%s/init+%ds.grb2", output_directory, (int) (t_write - t_init));
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "%s/init+%ds.grb2", output_directory, (int) (t_write - t_init));
    char *SAMPLE_FILENAME = "./grids/grib_template.grb2";
    FILE *SAMPLE_FILE;
    int err = 0;
    int retval;
    if (t_init < 0)
    	exit(1);
    FILE *OUT_GRIB;
    OUT_GRIB = fopen(OUTPUT_FILE, "w+");
    double *temperature = malloc(NO_OF_SCALARS*sizeof(double));
    double *pot_temperature = malloc(NO_OF_SCALARS*sizeof(double));
    // Grib requires everything to be on horizontal levels.
    double *temperature_h = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *pressure_h = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *rh_h = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *wind_u_h = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *wind_v_h = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *rel_vort_h = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *wind_w_h = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *mslp = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *surface_p = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *t2 = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *tcdc = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *rprate = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *sprate = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *cape = malloc(NO_OF_SCALARS_H*sizeof(double));
    Curl_field *rel_vort = calloc(1, sizeof(Curl_field));
    double pressure_value, mslp_factor, surface_p_factor, temp_lowest_layer, temp_mslp, temp_surface, wind_0, wind_1, wind_u, wind_v, delta_z_temp, temp_gradient, temp_upper, temp_lower;
    double standard_vert_lapse_rate = 0.0065;
    long unsigned length = 4;
    codes_handle *handle_temperature_h = NULL;
    codes_handle *handle_pressure_h = NULL;
    codes_handle *handle_wind_u_h = NULL;
    codes_handle *handle_wind_v_h = NULL;
    codes_handle *handle_wind_w_h = NULL;
    codes_handle *handle_wind_u_10m_mean = NULL;
    codes_handle *handle_wind_v_10m_mean = NULL;
    codes_handle *handle_mslp = NULL;
    codes_handle *handle_surface_p = NULL;
    codes_handle *handle_t2 = NULL;
    codes_handle *handle_tcdc = NULL;
    codes_handle *handle_rprate = NULL;
    codes_handle *handle_sprate = NULL;
    codes_handle *handle_wind_10m_gusts = NULL;
    codes_handle *handle_rel_vort = NULL;
    codes_handle *handle_rh = NULL;
    codes_handle *handle_cape = NULL;
	temperature_diagnostics(state_write_out, temperature);
	pot_temp_diagnostics(temperature, state_write_out -> density_dry, state_write_out -> tracer_densities, pot_temperature);
	calc_rel_vort(state_write_out -> velocity_gas, *rel_vort, grid, dualgrid);
	int layer_index;
	double z_height, delta_z, cape_integrand, theta_prime, theta;
	double z_tropopause = 15e3;
    for (int i = 0; i < NO_OF_LAYERS; ++i)
    {
        for (int j = 0; j < NO_OF_SCALARS_H; ++j)
        {
        	temperature_h[j] = temperature[i*NO_OF_SCALARS_H + j];
        	pressure_h[j] = state_write_out -> density_dry[j + i*NO_OF_SCALARS_H]*R_D*temperature_h[j];
        	rh_h[j] = 100*rel_humidity(temperature_h[j], state_write_out -> tracer_densities[2*NO_OF_SCALARS + i*NO_OF_SCALARS_H + j]);
        }
        if (i == NO_OF_LAYERS - 1)
        {
            for (int j = 0; j < NO_OF_SCALARS_H; ++j)
            {
            	// Now the aim is to determine the value of the MSLP.
                temp_lowest_layer = temperature[i*NO_OF_SCALARS_H + j];
                pressure_value = state_write_out -> density_dry[j + i*NO_OF_SCALARS_H]*R_D*temp_lowest_layer;
                temp_mslp = temp_lowest_layer + standard_vert_lapse_rate*grid -> z_scalar[j + i*NO_OF_SCALARS_H];
                mslp_factor = pow(1 - (temp_mslp - temp_lowest_layer)/temp_mslp, grid -> gravity_m[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + j]/(R_D*standard_vert_lapse_rate));
                mslp[j] = pressure_value/mslp_factor;
				// Now the aim is to determine the value of the surface pressure.
				temp_surface = temp_lowest_layer + standard_vert_lapse_rate*(grid -> z_scalar[j + i*NO_OF_SCALARS_H] - grid -> z_vector[NO_OF_VECTORS - NO_OF_VECTORS_V + j]);
                surface_p_factor = pow(1 - (temp_surface - temp_lowest_layer)/temp_surface, grid -> gravity_m[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + j]/(R_D*standard_vert_lapse_rate));
				surface_p[j] = pressure_value/surface_p_factor;
				// Now the aim is to calculate the 2 m temperature.
                delta_z_temp = 2 - grid -> z_scalar[j + i*NO_OF_SCALARS_H];
                temp_upper = temperature[(i - 1)*NO_OF_SCALARS_H + j];
                temp_lower = temp_lowest_layer;
                temp_gradient = (temp_upper - temp_lower)/(grid -> z_scalar[j + (i - 1)*NO_OF_SCALARS_H] - grid -> z_scalar[j + i*NO_OF_SCALARS_H]);
                // Finally the temperature in 2 m height AGL can bo obtained via linear extrapolation.
                t2[j] = temp_lowest_layer + delta_z_temp*temp_gradient;
                z_height = grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + j];
                cape[j] = 0;
                theta_prime = pot_temperature[i*NO_OF_SCALARS_H + j];
                layer_index = i;
                cape[j] = 0;
                while (z_height < z_tropopause)
                {
	                theta = pot_temperature[layer_index*NO_OF_SCALARS_H + j];
                	delta_z = grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + j] - grid -> z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + j];
                	z_height += delta_z;
                	cape_integrand = grid -> gravity_m[i*NO_OF_VECTORS_PER_LAYER + j]*(theta_prime - theta)/theta;
                	if (cape_integrand > 0)
                		cape[j] += cape_integrand*delta_z;
                	--layer_index;
                }
                // Now come the hydrometeors.
                sprate[j] = 0;
                rprate[j] = 0;
                tcdc[j] = 0;
                for (int k = 0; k < NO_OF_CONDENSATED_TRACERS; ++k)
                {
                    for (int l = 0; l < NO_OF_LAYERS; ++l)
                    {
                        if (state_write_out -> tracer_densities[k*NO_OF_SCALARS + l*NO_OF_SCALARS_H + j] > 0)
                            tcdc[j] = 100*1;
                    }
                }
                for (int k = 0; k < NO_OF_SOLID_TRACERS; ++k)
                    sprate[j] += fmax(ret_sink_velocity(k, 0, 0.001)*state_write_out -> tracer_densities[k*NO_OF_SCALARS + i*NO_OF_SCALARS_H + j], 0);
                for (int k = NO_OF_SOLID_TRACERS; k < NO_OF_CONDENSATED_TRACERS; ++k)
                    rprate[j] += fmax(ret_sink_velocity(k, 0, 0.001)*state_write_out -> tracer_densities[k*NO_OF_SCALARS + i*NO_OF_SCALARS_H + j], 0);
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
            if ((retval = codes_set_double_array(handle_mslp, "values", mslp, NO_OF_SCALARS_H)))
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
            if ((retval = codes_set_long(handle_surface_p, "typeOfFirstFixedSurface", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "scaledValueOfFirstFixedSurface", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "scaleFactorOfFirstFixedSurface", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_surface_p, "level", 0)))
                ECCERR(retval);
            if ((retval = codes_set_double_array(handle_surface_p, "values", surface_p, NO_OF_SCALARS_H)))
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
            if ((retval = codes_set_double_array(handle_t2, "values", t2, NO_OF_SCALARS_H)))
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
            if ((retval = codes_set_double_array(handle_rprate, "values", rprate, NO_OF_SCALARS_H)))
                ECCERR(retval);
            if ((retval = codes_write_message(handle_rprate, OUTPUT_FILE, "a")))
                ECCERR(retval);
			codes_handle_delete(handle_rprate);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_cape = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
            if ((retval = codes_set_long(handle_cape, "discipline", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "centre", 255)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "significanceOfReferenceTime", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "productionStatusOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "typeOfProcessedData", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "indicatorOfUnitOfTimeRange", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "stepUnits", 13)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "dataDate", data_date)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "dataTime", data_time)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "forecastTime", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "stepRange", t_write - t_init)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "typeOfGeneratingProcess", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "parameterCategory", 7)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "parameterNumber", 6)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "typeOfFirstFixedSurface", 103)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "scaledValueOfFirstFixedSurface", 0)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "scaleFactorOfFirstFixedSurface", 1)))
                ECCERR(retval);
            if ((retval = codes_set_long(handle_cape, "level", 0)))
                ECCERR(retval);
            if ((retval = codes_set_double_array(handle_cape, "values", cape, NO_OF_SCALARS_H)))
                ECCERR(retval);
            if ((retval = codes_write_message(handle_cape, OUTPUT_FILE, "a")))
                ECCERR(retval);
			codes_handle_delete(handle_cape);
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
            if ((retval = codes_set_double_array(handle_sprate, "values", sprate, NO_OF_SCALARS_H)))
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
            if ((retval = codes_set_double_array(handle_tcdc, "values", tcdc, NO_OF_SCALARS_H)))
                ECCERR(retval);
            if ((retval = codes_write_message(handle_tcdc, OUTPUT_FILE, "a")))
                ECCERR(retval);
			codes_handle_delete(handle_tcdc);
        }
    	SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_temperature_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
		    ECCERR(err);
		fclose(SAMPLE_FILE);
        if ((retval = codes_set_long(handle_temperature_h, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "stepUnits", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "dataDate", data_date)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "dataTime", data_time)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "forecastTime", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "stepRange", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "parameterCategory", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "parameterNumber", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "typeOfFirstFixedSurface", 26)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "scaledValueOfFirstFixedSurface", i)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_temperature_h, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_temperature_h, "values", temperature_h, NO_OF_SCALARS_H)))
            ECCERR(retval);
        if (i == 0)
        {
            if ((retval = codes_write_message(handle_temperature_h, OUTPUT_FILE, "w")))
                ECCERR(retval);
        }
        else
        {
            if ((retval = codes_write_message(handle_temperature_h, OUTPUT_FILE, "a")))
                ECCERR(retval);
        }
		codes_handle_delete(handle_temperature_h);
		SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_pressure_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
		    ECCERR(err);
		fclose(SAMPLE_FILE);
        if ((retval = codes_set_long(handle_pressure_h, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "stepUnits", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "dataDate", data_date)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "dataTime", data_time)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "forecastTime", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "stepRange", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "parameterCategory", 3)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "parameterNumber", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "typeOfFirstFixedSurface", 26)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "scaledValueOfFirstFixedSurface", i)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pressure_h, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_pressure_h, "values", pressure_h, NO_OF_SCALARS_H)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_pressure_h, OUTPUT_FILE, "a")))
            ECCERR(retval);
		codes_handle_delete(handle_pressure_h);
		SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_rh = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
		    ECCERR(err);
		fclose(SAMPLE_FILE);
        if ((retval = codes_set_long(handle_rh, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "stepUnits", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "dataDate", data_date)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "dataTime", data_time)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "forecastTime", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "stepRange", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "parameterCategory", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "parameterNumber", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "typeOfFirstFixedSurface", 26)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "scaledValueOfFirstFixedSurface", i)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rh, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_rh, "values", rh_h, NO_OF_SCALARS_H)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_rh, OUTPUT_FILE, "a")))
            ECCERR(retval);
		codes_handle_delete(handle_rh);
        for (int j = 0; j < NO_OF_VECTORS_H; ++j)
        {
            wind_0 = state_write_out -> velocity_gas[j + i*NO_OF_VECTORS_H + (i + 1)*NO_OF_VECTORS_V];
            recov_hor_par_pri(state_write_out -> velocity_gas, i, j, &wind_1, grid);
            passive_turn(wind_0, wind_1, -grid -> direction[j], &wind_u, &wind_v);
            wind_u_h[j] = wind_u;
            wind_v_h[j] = wind_v;
            rel_vort_h[j] = (*rel_vort)[NO_OF_DUAL_VECTORS_H + i*(NO_OF_VECTORS_H + NO_OF_DUAL_VECTORS_H) + j];
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
        if ((retval = codes_set_double_array(handle_wind_u_h, "values", wind_u_h, NO_OF_VECTORS_H)))
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
        if ((retval = codes_set_double_array(handle_wind_v_h, "values", wind_v_h, NO_OF_VECTORS_H)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_wind_v_h, OUTPUT_FILE, "a")))
            ECCERR(retval);
		codes_handle_delete(handle_wind_v_h);
		SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_rel_vort = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
		    ECCERR(err);
		fclose(SAMPLE_FILE);
        if ((retval = codes_set_long(handle_rel_vort, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "stepUnits", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "dataDate", data_date)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "dataTime", data_time)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "forecastTime", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "stepRange", t_write - t_init)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "parameterCategory", 2)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "parameterNumber", 12)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "typeOfFirstFixedSurface", 26)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "scaledValueOfFirstFixedSurface", i)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "scaleFactorOfFirstFixedSurface", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_rel_vort, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_rel_vort, "values", rel_vort_h, NO_OF_VECTORS_H)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_rel_vort, OUTPUT_FILE, "a")))
            ECCERR(retval);
		codes_handle_delete(handle_rel_vort);
    }
	double *wind_10_m_both_dir_array = malloc(2*min_no_of_output_steps*NO_OF_VECTORS_H*sizeof(double));
	double wind_10_m_downscale_factor = 0.667;
	double wind_tangential;
	int time_step_10_m_wind, h_index;
	double *wind_10_m_speed = malloc(min_no_of_output_steps*NO_OF_VECTORS_H*sizeof(double));
	double *wind_10_m_mean_u = malloc(NO_OF_VECTORS_H*sizeof(double));
	double *wind_10_m_mean_v = malloc(NO_OF_VECTORS_H*sizeof(double));
	for (int j = 0; j < min_no_of_output_steps*NO_OF_VECTORS_H; ++j)
	{
		time_step_10_m_wind = j/NO_OF_VECTORS_H;
		h_index = j - time_step_10_m_wind*NO_OF_VECTORS_H;
	    wind_10_m_both_dir_array[2*j + 0] = wind_10_m_downscale_factor*wind_h_lowest_layer_array[j];
	    wind_tangential = 0;
		for (int i = 0; i < 10; ++i)
			wind_tangential += grid -> trsk_modified_weights[10*h_index + i]*wind_h_lowest_layer_array[time_step_10_m_wind*NO_OF_VECTORS_H + grid -> trsk_modified_velocity_indices[10*h_index + i]];
	    wind_10_m_both_dir_array[2*j + 1] = wind_10_m_downscale_factor*wind_tangential;
	    wind_10_m_speed[j] = sqrt(pow(wind_10_m_both_dir_array[2*j + 0], 2) + pow(wind_10_m_both_dir_array[2*j + 1], 2));
	    if (time_step_10_m_wind == 0)
	    {
	    	wind_10_m_mean_u[h_index] = 1.0/min_no_of_output_steps*wind_10_m_both_dir_array[2*j + 0];
	    	wind_10_m_mean_v[h_index] = 1.0/min_no_of_output_steps*wind_10_m_both_dir_array[2*j + 1];
	    }
	    else
	    {
	    	wind_10_m_mean_u[h_index] += 1.0/min_no_of_output_steps*wind_10_m_both_dir_array[2*j + 0];
	    	wind_10_m_mean_v[h_index] += 1.0/min_no_of_output_steps*wind_10_m_both_dir_array[2*j + 1];
	    }
	}
	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
	{
        passive_turn(wind_10_m_mean_u[i], wind_10_m_mean_v[i], -grid -> direction[i], &wind_u, &wind_v);
        wind_10_m_mean_u[i] = wind_u;
        wind_10_m_mean_v[i] = wind_v;
	}
	double standard_deviation;
	double gusts_parameter = 3;
	double *wind_10_m_gusts_speed = malloc(NO_OF_VECTORS_H*sizeof(double));
	double *vector_for_std_deviation = malloc(min_no_of_output_steps*sizeof(double));
	double wind_speed_10_m_mean;
	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		wind_speed_10_m_mean = 0;
		for (int j = 0; j < min_no_of_output_steps; ++j)
		{
			vector_for_std_deviation[j] = wind_10_m_speed[j*NO_OF_VECTORS_H + i];
			wind_speed_10_m_mean += 1.0/min_no_of_output_steps*wind_10_m_speed[j*NO_OF_VECTORS_H + i];
		}
		standard_deviation = calc_std_dev(vector_for_std_deviation, min_no_of_output_steps);
		if (t_write != t_init)
			wind_10_m_gusts_speed[i] = wind_speed_10_m_mean + gusts_parameter*standard_deviation;
		else
			wind_10_m_gusts_speed[i] = (1 + 0.2)*wind_speed_10_m_mean;
	}
	free(vector_for_std_deviation);
	SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
	handle_wind_u_10m_mean = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	fclose(SAMPLE_FILE);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "discipline", 0)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "centre", 255)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "significanceOfReferenceTime", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "productionStatusOfProcessedData", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "typeOfProcessedData", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "indicatorOfUnitOfTimeRange", 13)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "stepUnits", 13)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "dataDate", data_date)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "dataTime", data_time)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "forecastTime", t_write - t_init)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "stepRange", t_write - t_init)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "typeOfGeneratingProcess", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "parameterCategory", 2)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "parameterNumber", 2)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "typeOfFirstFixedSurface", 103)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "scaledValueOfFirstFixedSurface", 10)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "scaleFactorOfFirstFixedSurface", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_u_10m_mean, "level", 10)))
	    ECCERR(retval);
	if ((retval = codes_set_double_array(handle_wind_u_10m_mean, "values", wind_10_m_mean_u, NO_OF_VECTORS_H)))
	    ECCERR(retval);
	if ((retval = codes_write_message(handle_wind_u_10m_mean, OUTPUT_FILE, "a")))
	    ECCERR(retval);
	codes_handle_delete(handle_wind_u_10m_mean);
	free(wind_10_m_mean_u);
	SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
	handle_wind_v_10m_mean = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	fclose(SAMPLE_FILE);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "discipline", 0)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "centre", 255)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "significanceOfReferenceTime", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "productionStatusOfProcessedData", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "typeOfProcessedData", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "indicatorOfUnitOfTimeRange", 13)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "stepUnits", 13)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "dataDate", data_date)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "dataTime", data_time)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "forecastTime", t_write - t_init)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "stepRange", t_write - t_init)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "typeOfGeneratingProcess", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "parameterCategory", 2)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "parameterNumber", 3)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "typeOfFirstFixedSurface", 103)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "scaledValueOfFirstFixedSurface", 10)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "scaleFactorOfFirstFixedSurface", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_v_10m_mean, "level", 10)))
	    ECCERR(retval);
	if ((retval = codes_set_double_array(handle_wind_v_10m_mean, "values", wind_10_m_mean_v, NO_OF_VECTORS_H)))
	    ECCERR(retval);
	if ((retval = codes_write_message(handle_wind_v_10m_mean, OUTPUT_FILE, "a")))
	    ECCERR(retval);
	free(wind_10_m_mean_v);
	codes_handle_delete(handle_wind_v_10m_mean);
	SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
	handle_wind_10m_gusts = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	fclose(SAMPLE_FILE);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "discipline", 0)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "centre", 255)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "significanceOfReferenceTime", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "productionStatusOfProcessedData", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "typeOfProcessedData", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "indicatorOfUnitOfTimeRange", 13)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "stepUnits", 13)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "dataDate", data_date)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "dataTime", data_time)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "forecastTime", t_write - t_init)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "stepRange", t_write - t_init)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "typeOfGeneratingProcess", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "parameterCategory", 2)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "parameterNumber", 22)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "typeOfFirstFixedSurface", 103)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "scaledValueOfFirstFixedSurface", 10)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "scaleFactorOfFirstFixedSurface", 1)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle_wind_10m_gusts, "level", 10)))
	    ECCERR(retval);
	if ((retval = codes_set_double_array(handle_wind_10m_gusts, "values", wind_10_m_gusts_speed, NO_OF_VECTORS_H)))
	    ECCERR(retval);
	if ((retval = codes_write_message(handle_wind_10m_gusts, OUTPUT_FILE, "a")))
	    ECCERR(retval);
	codes_handle_delete(handle_wind_10m_gusts);
	free(wind_10_m_speed);
	free(wind_10_m_gusts_speed);
	free(wind_10_m_both_dir_array);
    free(t2);
    free(mslp);
    free(surface_p);
    free(rprate);
    free(sprate);
    free(tcdc);
    free(pot_temperature);
    free(temperature);
    free(rel_vort_h);
    free(rel_vort);
    free(temperature_h);
    free(pressure_h);
    free(wind_u_h);
    free(wind_v_h);
    free(rh_h);
    free(cape);
	SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
	handle_wind_w_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
	    ECCERR(err);
	fclose(SAMPLE_FILE);
    for (int i = 0; i < NO_OF_LEVELS; ++i)
    {
        for (int j = 0; j < NO_OF_VECTORS_V; j++)
            wind_w_h[j] = state_write_out -> velocity_gas[j + i*NO_OF_VECTORS_PER_LAYER];
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
        if ((retval = codes_set_double_array(handle_wind_w_h, "values", wind_w_h, NO_OF_VECTORS_V)))
            ECCERR(retval);
        if ((retval = codes_write_message(handle_wind_w_h, OUTPUT_FILE, "a")))
            ECCERR(retval);
    }
    codes_handle_delete(handle_wind_w_h);
    free(wind_w_h);
    fclose(OUT_GRIB);
    return 0;
}

int write_out_integral(State *state_write_out, double t_write, char output_directory[], Grid *grid, Dualgrid *dualgrid, int integral_id)
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
   		sprintf(INTEGRAL_FILE_PRE, "%s/%s", output_directory, "dry_mass");
    if (integral_id == 1)
   		sprintf(INTEGRAL_FILE_PRE, "%s/%s", output_directory, "entropy");
    if (integral_id == 2)
   		sprintf(INTEGRAL_FILE_PRE, "%s/%s", output_directory, "energy");
    INTEGRAL_FILE_LENGTH = strlen(INTEGRAL_FILE_PRE);
    char *INTEGRAL_FILE = malloc((INTEGRAL_FILE_LENGTH + 1)*sizeof(char));
    sprintf(INTEGRAL_FILE, "%s", INTEGRAL_FILE_PRE);
    free(INTEGRAL_FILE_PRE);
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
    	kinetic_energy(state_write_out -> velocity_gas, *e_kin_density, grid);
    	scalar_times_scalar(state_write_out -> density_dry, *e_kin_density, *e_kin_density);
    	global_scalar_integrator(*e_kin_density, grid, &kinetic_integral);
    	free(e_kin_density);
    	Scalar_field *pot_energy_density = malloc(sizeof(Scalar_field));
    	scalar_times_scalar(state_write_out -> density_dry, grid -> gravity_potential, *pot_energy_density);
    	global_scalar_integrator(*pot_energy_density, grid, &potential_integral);
    	free(pot_energy_density);
    	Scalar_field *int_energy_density = malloc(sizeof(Scalar_field));
    	temperature_diagnostics(state_write_out, *int_energy_density);
    	scalar_times_scalar(state_write_out -> density_dry, *int_energy_density, *int_energy_density);
    	global_scalar_integrator(*int_energy_density, grid, &internal_integral);
    	fprintf(global_integral_file, "%lf\t%lf\t%lf\t%lf\n", t_write, kinetic_integral, potential_integral, C_D_V*internal_integral);
    	free(int_energy_density);
    	fclose(global_integral_file);
    }
    free(INTEGRAL_FILE);
	return 0;
}

double calc_std_dev(double vector_for_std_deviation[], int no_of_values)
{
	double mean = 0;
	for (int i = 0; i < no_of_values; ++i)
		mean += 1.0/no_of_values*vector_for_std_deviation[i];
	double result = 0;
	for (int i = 0; i < no_of_values; ++i)
		result += pow(vector_for_std_deviation[i] - mean, 2);
	result = 1/sqrt(no_of_values)*sqrt(result);
	return result;
}

