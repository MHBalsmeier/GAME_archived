/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
Test states can be generated with this code. TEST_ID table:
0:	standard atmosphere
1:	standard atmosphere with Gaussian mountain
2:	JW dry unperturbed
3:	JW dry perturbed
4:	JW moist unperturbed
5:	JW moist perturbed
For more specific details see handbook.
*/

#include <stdlib.h>
#include "enum.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include "eccodes.h"
#include "geos95.h"
#include "atmostracers.h"
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define P_0 100000
#define OMEGA (7.292115e-5)
#define C_D_P 1005.0

// constants specifying the grid
const double TOA = 30000;
const double G = 9.80616;

// constants needed for the JW test state
const double TROPO_HEIGHT = 12e3;
const double T_SFC = 273.15 + 15;
double T_0 = 288;
const double TEMP_GRADIENT = -0.65/100;
double GAMMA = 0.005;
const double DELTA_T = 4.8e5;
const double ETA_T = 0.2;
const double U_0 = 35;
const double ETA_0 = 0.252;

int find_pressure_value(double, double, double *);
int find_z_from_p(double, double, double *);

int main(int argc, char *argv[])
{
	int TEST_ID;
   	TEST_ID = strtod(argv[1], NULL);
   	// determining the orography ID as a function of the test ID
	int ORO_ID;
	if (TEST_ID == 0)
		ORO_ID = 0;
	if (TEST_ID == 1)
		ORO_ID = 1;
	if (TEST_ID > 1)
		ORO_ID = 2;
    double *direction = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *latitude_scalar = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *latitude_vector = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *longitude_vector = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *z_scalar = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *z_vector = malloc(NUMBER_OF_VECTORS*sizeof(double));
    double *gravity_potential = malloc(NUMBER_OF_VECTORS*sizeof(double));
    int ncid, retval;
    short GEO_PROP_FILE_LENGTH = 100;
    char *GEO_PROP_FILE_PRE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE_PRE, "../grid_generator/nc_files/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NUMBER_OF_LAYERS, (int) TOA, ORO_ID, NUMBER_OF_ORO_LAYERS);
    GEO_PROP_FILE_LENGTH = strlen(GEO_PROP_FILE_PRE);
    free(GEO_PROP_FILE_PRE);
    char *GEO_PROP_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE, "../grid_generator/nc_files/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NUMBER_OF_LAYERS, (int) TOA, ORO_ID, NUMBER_OF_ORO_LAYERS);
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid)))
        NCERR(retval);
    free(GEO_PROP_FILE);
    int direction_id, latitude_scalar_id, longitude_scalar_id, latitude_vector_id, longitude_vector_id, z_scalar_id, z_vector_id, gravity_potential_id;
    if ((retval = nc_inq_varid(ncid, "direction", &direction_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "latitude_scalar", &latitude_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_scalar", &longitude_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "latitude_vector", &latitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_vector", &longitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_scalar", &z_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_vector", &z_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "gravity_potential", &gravity_potential_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, direction_id, &direction[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_scalar_id, &latitude_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_scalar_id, &longitude_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_vector_id, &latitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_vector_id, &longitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, z_scalar_id, &z_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, z_vector_id, &z_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, gravity_potential_id, &gravity_potential[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
        NCERR(retval);
    char *SAMPLE_FILE = "grib_files/grib_template.grb2";
    FILE *SAMPLE;
    int err = 0;
    short OUTPUT_FILE_LENGTH = 100;
    char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE_PRE, "grib_files/test_%d_B%dL%dT%d_O%d_OL%d_SCVT.grb2", TEST_ID, RES_ID, NUMBER_OF_LAYERS, (int) TOA, ORO_ID, NUMBER_OF_ORO_LAYERS);
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "grib_files/test_%d_B%dL%dT%d_O%d_OL%d_SCVT.grb2", TEST_ID, RES_ID, NUMBER_OF_LAYERS, (int) TOA, ORO_ID, NUMBER_OF_ORO_LAYERS);
    codes_handle *handle_pot_temperature = NULL;
    codes_handle *handle_density = NULL;
    codes_handle *handle_wind_h = NULL;
    codes_handle *handle_wind_v = NULL;
    codes_handle *handle_water_vapour_density = NULL;
    codes_handle *handle_liquid_water_density = NULL;
    codes_handle *handle_solid_water_density = NULL;
    codes_handle *handle_liquid_water_temp = NULL;
    codes_handle *handle_solid_water_temp = NULL;
    // 3D arrays
    double *pressure = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *pot_temperature = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *temperature = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *rho = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *rel_humidity = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *water_vapour_density = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *liquid_water_density = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *solid_water_density = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *liquid_water_temp = malloc(NUMBER_OF_SCALARS*sizeof(double));
    double *solid_water_temp = malloc(NUMBER_OF_SCALARS*sizeof(double));
    // grib requires everything on horizontal levels
    double *pressure_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *pot_temperature_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *temperature_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *rho_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *water_vapour_density_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *liquid_water_density_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *solid_water_density_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *liquid_water_temp_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *solid_water_temp_h = malloc(NUMBER_OF_SCALARS_H*sizeof(double));
    double *wind_h = malloc(NUMBER_OF_VECTORS_H*sizeof(double));
    double *wind_v = malloc(NUMBER_OF_VECTORS_V*sizeof(double));
    const double TROPO_TEMP = T_SFC + TROPO_HEIGHT*TEMP_GRADIENT;
    double z_height;
    SAMPLE = fopen(SAMPLE_FILE, "r");
    handle_pot_temperature = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE);
    SAMPLE = fopen(SAMPLE_FILE, "r");
    handle_density = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE);
    SAMPLE = fopen(SAMPLE_FILE, "r");
    handle_water_vapour_density = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE);
    SAMPLE = fopen(SAMPLE_FILE, "r");
    handle_liquid_water_density = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE);
    SAMPLE = fopen(SAMPLE_FILE, "r");
    handle_solid_water_density = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE);
    SAMPLE = fopen(SAMPLE_FILE, "r");
    handle_liquid_water_temp = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE);
    SAMPLE = fopen(SAMPLE_FILE, "r");
    handle_solid_water_temp = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE);
    SAMPLE = fopen(SAMPLE_FILE, "r");
    handle_wind_h = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE);
    double lat, lon;
    double u, eta, eta_v, T_perturb, distance, pressure_value;
    double u_p = 1.0;
    double distance_scale = RADIUS/10;
    double lat_perturb = 2*M_PI/9;
    double lon_perturb = M_PI/9;
    int layer_index, h_index;
    // 3D scalar fields determined here, apart from density
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
    	layer_index = i/NUMBER_OF_SCALARS_H;
    	h_index = i - layer_index*NUMBER_OF_SCALARS_H;
        lat = latitude_scalar[h_index];
        lon = longitude_scalar[h_index];
        z_height = z_scalar[i];
        rel_humidity[i] = 0;
        // standard atmosphere
        if (TEST_ID == 0 || TEST_ID == 1)
        {
            if (z_height < TROPO_HEIGHT)
            {
                temperature[i] = T_SFC + z_height*TEMP_GRADIENT;
                pressure[i] = 101325*pow(1 + TEMP_GRADIENT*z_height/T_SFC, -G/(R_D*TEMP_GRADIENT));
            }
            else
            {
                temperature[i] = TROPO_TEMP;
                pressure[i] = 101325*pow(1 + TEMP_GRADIENT*TROPO_HEIGHT/T_SFC, -G/(R_D*TEMP_GRADIENT))*exp(-G*(z_height - TROPO_HEIGHT)/(R_D*TROPO_TEMP));
            }
        }
        // JW atmosphere
        if (TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5)
        {
            find_pressure_value(lat, z_height, &pressure_value);
            pressure[i] = pressure_value;
            eta = pressure[i]/P_0;
            eta_v = (eta - ETA_0)*M_PI/2;
            T_perturb = 3.0/4.0*eta*M_PI*U_0/R_D*sin(eta_v)*pow(cos(eta_v), 0.5)*((-2*pow(sin(lat), 6)*(pow(cos(lat), 2) + 1.0/3.0) + 10.0/63.0)*2*U_0*pow(cos(eta_v), 1.5) + RADIUS*OMEGA*(8.0/5.0*pow(cos(lat), 3)*(pow(sin(lat), 2) + 2.0/3.0) - M_PI/4.0));
            if (eta >= ETA_T)
            {
                temperature[i] = T_0*pow(eta, R_D*GAMMA/G) + T_perturb;
                if (TEST_ID == 4 || TEST_ID == 5)
                    rel_humidity[i] = 0.7;
                else
                    rel_humidity[i] = 0;
            }
            else
                temperature[i] = T_0*pow(eta, R_D*GAMMA/G) + DELTA_T*pow(ETA_T - eta, 5) + T_perturb;
        }
        liquid_water_density[i] = 0;
        solid_water_density[i] = 0;
        liquid_water_temp[i] = temperature[i];
        solid_water_temp[i] = temperature[i];
    }
    // density is determined out of the hydrostatic equation
    double entropy_value, temperature_mean, delta_temperature, delta_gravity_potential, pot_temp_value, lower_entropy_value;
    for (int i = NUMBER_OF_SCALARS - 1; i >= 0; --i)
    {
    	layer_index = i/NUMBER_OF_SCALARS_H;
    	// at the lowest layer the density is set using the equation of state, can be considered a boundary condition
    	if (layer_index == NUMBER_OF_LAYERS - 1)
    	{
        	pot_temperature[i] = temperature[i]*pow(P_0/pressure[i], R_D/C_D_P);
        	rho[i] = pressure[i]/(R_D*temperature[i]);
        }
        else
        {
        	lower_entropy_value = C_D_P*log(temperature[i + NUMBER_OF_SCALARS_H]*pow(P_0/(rho[i + NUMBER_OF_SCALARS_H]*R_D*temperature[i + NUMBER_OF_SCALARS_H]), R_D/C_D_P));
        	temperature_mean = 0.5*(temperature[i] + temperature[i + NUMBER_OF_SCALARS_H]);
        	delta_temperature = temperature[i] - temperature[i + NUMBER_OF_SCALARS_H];
        	delta_gravity_potential = gravity_potential[i] - gravity_potential[i + NUMBER_OF_SCALARS_H];
        	entropy_value = lower_entropy_value + (delta_gravity_potential + C_D_P*delta_temperature)/temperature_mean;
        	pot_temp_value = exp(entropy_value/C_D_P);
        	pot_temperature[i] = pot_temp_value;
        	pressure_value = P_0*pow(temperature[i]/pot_temp_value, C_D_P/R_D);
        	rho[i] = pressure_value/(R_D*temperature[i]);
        }
        water_vapour_density[i] = water_vapour_density_from_rel_humidity(rel_humidity[i], temperature[i], rho[i]);
        if (water_vapour_density[i] < 0)
        	printf("water_vapour_density negative.\n.");
    }
    free(gravity_potential);
    for (int i = 0; i < NUMBER_OF_LAYERS; ++i)
    {
    	// Grib wants everything as 2D levels. These arrays are constructed here out of the 3D arrays.
        for (int j = 0; j < NUMBER_OF_SCALARS_H; ++j)
        {
            rho_h[j] = rho[i*NUMBER_OF_SCALARS_H + j];
            pot_temperature_h[j] = pot_temperature[i*NUMBER_OF_SCALARS_H + j];
            water_vapour_density_h[j] = water_vapour_density[i*NUMBER_OF_SCALARS_H + j];
            liquid_water_density_h[j] = liquid_water_density[i*NUMBER_OF_SCALARS_H + j];
            solid_water_density_h[j] = solid_water_density[i*NUMBER_OF_SCALARS_H + j];
            liquid_water_temp_h[j] = liquid_water_temp[i*NUMBER_OF_SCALARS_H + j];
            solid_water_temp_h[j] = solid_water_temp[i*NUMBER_OF_SCALARS_H + j];
        }
        if ((retval = codes_set_long(handle_pot_temperature, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "dataDate", 20000101)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "dataTime", 0000)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "parameterCategory", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "parameterNumber", 2)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "scaledValueOfFirstFixedSurface", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "scaleFactorOfFirstFixedSurface", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_pot_temperature, "level", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_pot_temperature, "values", pot_temperature_h, NUMBER_OF_SCALARS_H)))
            ECCERR(retval);
        if (i == 0)
            codes_write_message(handle_pot_temperature, OUTPUT_FILE, "w");
        else
            codes_write_message(handle_pot_temperature, OUTPUT_FILE, "a");
        if ((retval = codes_set_long(handle_density, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "dataDate", 20000101)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "dataTime", 0000)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "parameterCategory", 3)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "parameterNumber", 10)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "scaledValueOfFirstFixedSurface", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "scaleFactorOfFirstFixedSurface", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_density, "level", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_density, "values", rho_h, NUMBER_OF_SCALARS_H)))
            ECCERR(retval);
        codes_write_message(handle_density, OUTPUT_FILE, "a");
        if ((retval = codes_set_long(handle_water_vapour_density, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "dataDate", 20000101)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "dataTime", 0000)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "parameterCategory", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "parameterNumber", 18)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "scaledValueOfFirstFixedSurface", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "scaleFactorOfFirstFixedSurface", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_water_vapour_density, "level", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_water_vapour_density, "values", water_vapour_density_h, NUMBER_OF_SCALARS_H)))
            ECCERR(retval);
        codes_write_message(handle_water_vapour_density, OUTPUT_FILE, "a");
        if ((retval = codes_set_long(handle_liquid_water_density, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "dataDate", 20000101)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "dataTime", 0000)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "parameterCategory", 6)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "parameterNumber", 38)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "scaledValueOfFirstFixedSurface", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "scaleFactorOfFirstFixedSurface", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_density, "level", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_liquid_water_density, "values", liquid_water_density_h, NUMBER_OF_SCALARS_H)))
            ECCERR(retval);
        codes_write_message(handle_liquid_water_density, OUTPUT_FILE, "a");
        if ((retval = codes_set_long(handle_solid_water_density, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "dataDate", 20000101)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "dataTime", 0000)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "parameterCategory", 6)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "parameterNumber", 39)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "scaledValueOfFirstFixedSurface", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "scaleFactorOfFirstFixedSurface", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_density, "level", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_solid_water_density, "values", solid_water_density_h, NUMBER_OF_SCALARS_H)))
            ECCERR(retval);
        codes_write_message(handle_solid_water_density, OUTPUT_FILE, "a");
        if ((retval = codes_set_long(handle_liquid_water_temp, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "dataDate", 20000101)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "dataTime", 0000)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "parameterCategory", 6)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "parameterNumber", 192)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "scaledValueOfFirstFixedSurface", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "scaleFactorOfFirstFixedSurface", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_liquid_water_temp, "level", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_liquid_water_temp, "values", liquid_water_temp_h, NUMBER_OF_SCALARS_H)))
            ECCERR(retval);
        codes_write_message(handle_liquid_water_temp, OUTPUT_FILE, "a");
        if ((retval = codes_set_long(handle_solid_water_temp, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "dataDate", 20000101)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "dataTime", 0000)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "parameterCategory", 6)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "parameterNumber", 193)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "scaledValueOfFirstFixedSurface", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "scaleFactorOfFirstFixedSurface", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_solid_water_temp, "level", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_solid_water_temp, "values", solid_water_temp_h, NUMBER_OF_SCALARS_H)))
            ECCERR(retval);
        codes_write_message(handle_solid_water_temp, OUTPUT_FILE, "a");
        // wind fields are determind here
        for (int j = 0; j < NUMBER_OF_VECTORS_H; ++j)
        {
            lat = latitude_vector[j];
            lon = longitude_vector[j];
            z_height = z_vector[NUMBER_OF_VECTORS_V + j + i*NUMBER_OF_VECTORS_PER_LAYER];
            // standard atmosphere: no wind
            if (TEST_ID == 0 || TEST_ID == 1)
                wind_h[j] = 0;
            // JW test: specific wind field
            if (TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5)
            {
                find_pressure_value(lat, z_height, &pressure_value);
                eta = pressure_value/P_0;
                eta_v = (eta - ETA_0)*M_PI/2; 
                u = U_0*pow(cos(eta_v), 1.5)*pow(sin(2*lat), 2);
                if (TEST_ID == 3 || TEST_ID == 5)
                {
                    distance = calculate_distance_h(lat, lon, lat_perturb, lon_perturb, RADIUS);
                    u += u_p*exp(-pow(distance/distance_scale, 2));
                }
                wind_h[j] = u*cos(direction[j]);
            }
        }
        if ((retval = codes_set_long(handle_wind_h, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "dataDate", 20000101)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "dataTime", 0000)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "parameterCategory", 2)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "parameterNumber", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "scaledValueOfFirstFixedSurface", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "scaleFactorOfFirstFixedSurface", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_h, "level", round(z_height))))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_wind_h, "values", wind_h, NUMBER_OF_VECTORS_H)))
            ECCERR(retval);
        codes_write_message(handle_wind_h, OUTPUT_FILE, "a");
    }
    free(direction);
    free(latitude_vector);
    free(longitude_vector);
    free(pressure);
    free(pot_temperature);
    free(temperature);
    free(rho);
    free(water_vapour_density);
    free(liquid_water_density);
    free(solid_water_density);
    free(liquid_water_temp);
    free(solid_water_temp);
    free(z_scalar);
    free(rel_humidity);
    free(pressure_h);
    free(pot_temperature_h);
    free(temperature_h);
    free(rho_h);
    free(water_vapour_density_h);
    free(liquid_water_density_h);
    free(solid_water_density_h);
    free(liquid_water_temp_h);
    free(solid_water_temp_h);
    free(wind_h);
    codes_handle_delete(handle_pot_temperature);
    codes_handle_delete(handle_density);
    codes_handle_delete(handle_wind_h);
    codes_handle_delete(handle_water_vapour_density);
    codes_handle_delete(handle_liquid_water_density);
    codes_handle_delete(handle_solid_water_density);
    codes_handle_delete(handle_liquid_water_temp);
    codes_handle_delete(handle_solid_water_temp);
    SAMPLE = fopen(SAMPLE_FILE, "r");
    handle_wind_v = codes_handle_new_from_file(NULL, SAMPLE, PRODUCT_GRIB, &err);
    if (err != 0)
        ECCERR(err);
    fclose(SAMPLE);
    for (int i = 0; i < NUMBER_OF_LEVELS; ++i)
    {
        for (int j = 0; j < NUMBER_OF_VECTORS_V; ++j)
        {
            lat = latitude_scalar[j];
            lon = longitude_scalar[j];
            z_height = z_vector[j + i*NUMBER_OF_VECTORS_PER_LAYER];
            if (TEST_ID == 0 || TEST_ID == 1 || TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5)
                wind_v[j] = 0;
        }
        if ((retval = codes_set_long(handle_wind_v, "discipline", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "centre", 255)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "significanceOfReferenceTime", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "productionStatusOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "typeOfProcessedData", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "indicatorOfUnitOfTimeRange", 13)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "dataDate", 20000101)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "dataTime", 0000)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "typeOfGeneratingProcess", 1)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "parameterCategory", 2)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "parameterNumber", 9)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "typeOfFirstFixedSurface", 102)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "scaledValueOfFirstFixedSurface", z_height)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "scaleFactorOfFirstFixedSurface", 0)))
            ECCERR(retval);
        if ((retval = codes_set_long(handle_wind_v, "level", i)))
            ECCERR(retval);
        if ((retval = codes_set_double_array(handle_wind_v, "values", wind_v, NUMBER_OF_SCALARS_H)))
            ECCERR(retval);
        codes_write_message(handle_wind_v, OUTPUT_FILE, "a");
    }
    free(z_vector);
    free(latitude_scalar);
    free(longitude_scalar);
    codes_handle_delete(handle_wind_v);
    free(OUTPUT_FILE);
    free(wind_v);
    return 0;
}

int find_pressure_value(double lat, double z_height, double *result)
{
	// this function finds the pressure at a given height (as a function of latitude) for the JW test by iterative calls to the function find_z_from_P
    double p = P_0/2;
    double precision = 0.0001;
    double z;
    double current_max = P_0;
    double current_min = 0;
    find_z_from_p(lat, p, &z);
    while (fabs(z - z_height) > precision)
    {
        if (z < z_height)
        {
            current_max = p;
            p = 0.5*(current_min + p);
        }
        else
        {
            current_min = p;
            p = 0.5*(current_max + p);
        }
        find_z_from_p(lat, p, &z);
    }
    *result = p;
    return 0;
}

int find_z_from_p(double lat, double p, double *result)
{
	// this function converts a preessure value into a geomtrical height (as a function of latitude) for the JW test
    double z;
    double eta = p/P_0;
    double phi, phi_bg, phi_perturb;
    double eta_v = (eta - ETA_0)*M_PI/2;
    if (eta >= ETA_T)
        phi_bg = T_0*G/GAMMA*(1 - pow(eta, R_D*GAMMA/G));
    else
        phi_bg = T_0*G/GAMMA*(1 - pow(eta, R_D*GAMMA/G)) - R_D*DELTA_T*((log(eta/ETA_T) + 137.0/60.0)*pow(ETA_T, 5) - 5*eta*pow(ETA_T, 4) + 5*pow(ETA_T, 3)*pow(eta, 2) - 10.0/3.0*pow(ETA_T, 2)*pow(eta, 3) + 5.0/4.0*ETA_T*pow(eta, 4) - 1.0/5.0*pow(eta, 5));
    phi_perturb = U_0*pow(cos(eta_v), 1.5)*((-2*pow(sin(lat), 6)*(pow(cos(lat), 2) + 1.0/3.0) + 10.0/63.0)*U_0*pow(cos(eta_v), 1.5) + RADIUS*OMEGA*(8.0/5.0*pow(cos(lat), 3)*(pow(sin(lat), 2) + 2.0/3.0) - M_PI/4.0));
    phi = phi_bg + phi_perturb;
    z = phi/G;
    *result = z;
    return 0;
}

















