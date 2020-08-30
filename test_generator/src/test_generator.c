/*
This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
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
#include "geos95.h"
#include "atmostracers.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define P_0 100000.0
#define OMEGA (7.292115e-5)
#define C_D_P 1005.0

// constants specifying the grid
const double TOA = 30000;

// constants needed for the JW test state
const double G = 9.80616;
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
	if (TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5)
		ORO_ID = 2;
	if (TEST_ID == 6 || TEST_ID == 7)
		ORO_ID = 3;
    double *direction = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *latitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *latitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *longitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *z_scalar = malloc(NO_OF_SCALARS*sizeof(double));
    double *z_vector = malloc(NO_OF_VECTORS*sizeof(double));
    double *gravity_potential = malloc(NO_OF_VECTORS*sizeof(double));
    int ncid_grid, retval;
    int GEO_PROP_FILE_LENGTH = 100;
    char *GEO_PROP_FILE_PRE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE_PRE, "../grid_generator/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    GEO_PROP_FILE_LENGTH = strlen(GEO_PROP_FILE_PRE);
    free(GEO_PROP_FILE_PRE);
    char *GEO_PROP_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE, "../grid_generator/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid_grid)))
        NCERR(retval);
    free(GEO_PROP_FILE);
    int direction_id, latitude_scalar_id, longitude_scalar_id, latitude_vector_id, longitude_vector_id, z_scalar_id, z_vector_id, gravity_potential_id;
    if ((retval = nc_inq_varid(ncid_grid, "direction", &direction_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "latitude_scalar", &latitude_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_scalar", &longitude_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "latitude_vector", &latitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_vector", &longitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_scalar", &z_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_vector", &z_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "gravity_potential", &gravity_potential_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, direction_id, &direction[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitude_scalar_id, &latitude_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitude_scalar_id, &longitude_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitude_vector_id, &latitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitude_vector_id, &longitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_scalar_id, &z_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_vector_id, &z_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, gravity_potential_id, &gravity_potential[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid_grid)))
        NCERR(retval);
    int OUTPUT_FILE_LENGTH = 100;
    char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE_PRE, "test_states/test_%d_B%dL%dT%d_O%d_OL%d_SCVT.nc", TEST_ID, RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "test_states/test_%d_B%dL%dT%d_O%d_OL%d_SCVT.nc", TEST_ID, RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    double *pressure = malloc(NO_OF_SCALARS*sizeof(double));
    double *temperature = malloc(NO_OF_SCALARS*sizeof(double));
    double *rho = malloc(NO_OF_SCALARS*sizeof(double));
    double *rel_humidity = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *water_vapour_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_temp = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_temp = malloc(NO_OF_SCALARS*sizeof(double));
    const double TROPO_TEMP = T_SFC + TROPO_HEIGHT*TEMP_GRADIENT;
    double z_height;
    double lat, lon;
    double u, eta, eta_v, T_perturb, distance, pressure_value;
    double u_p = 1.0;
    double distance_scale = RADIUS/10;
    double lat_perturb = 2*M_PI/9;
    double lon_perturb = M_PI/9;
    int layer_index, h_index;
    // 3D scalar fields determined here, apart from density
    #pragma omp parallel for private(layer_index, h_index, lat, lon, z_height, eta, eta_v, T_perturb, pressure_value)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
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
        if (TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5 || TEST_ID == 6)
        {
            find_pressure_value(lat, z_height, &pressure_value);
            pressure[i] = pressure_value;
            eta = pressure[i]/P_0;
            eta_v = (eta - ETA_0)*M_PI/2;
            T_perturb = 3.0/4.0*eta*M_PI*U_0/R_D*sin(eta_v)*pow(cos(eta_v), 0.5)*((-2*pow(sin(lat), 6)*(pow(cos(lat), 2) + 1.0/3.0) + 10.0/63.0)*2*U_0*pow(cos(eta_v), 1.5) + RADIUS*OMEGA*(8.0/5.0*pow(cos(lat), 3)*(pow(sin(lat), 2) + 2.0/3.0) - M_PI/4.0));
            if (eta >= ETA_T)
            {
                temperature[i] = T_0*pow(eta, R_D*GAMMA/G) + T_perturb;
                if (TEST_ID == 4 || TEST_ID == 5 || TEST_ID == 7)
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
    for (int i = NO_OF_SCALARS - 1; i >= 0; --i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	// at the lowest layer the density is set using the equation of state, can be considered a boundary condition
    	if (layer_index == NO_OF_LAYERS - 1)
    	{
        	rho[i] = pressure[i]/(R_D*temperature[i]);
        }
        else
        {
        	lower_entropy_value = C_D_P*log(temperature[i + NO_OF_SCALARS_H]*pow(P_0/(rho[i + NO_OF_SCALARS_H]*R_D*temperature[i + NO_OF_SCALARS_H]), R_D/C_D_P));
        	temperature_mean = 0.5*(temperature[i] + temperature[i + NO_OF_SCALARS_H]);
        	delta_temperature = temperature[i] - temperature[i + NO_OF_SCALARS_H];
        	delta_gravity_potential = gravity_potential[i] - gravity_potential[i + NO_OF_SCALARS_H];
        	entropy_value = lower_entropy_value + (delta_gravity_potential + C_D_P*delta_temperature)/temperature_mean;
        	pot_temp_value = exp(entropy_value/C_D_P);
        	pressure_value = P_0*pow(temperature[i]/pot_temp_value, C_D_P/R_D);
        	rho[i] = pressure_value/(R_D*temperature[i]);
        }
        water_vapour_density[i] = water_vapour_density_from_rel_humidity(rel_humidity[i], temperature[i], rho[i]);
        if (water_vapour_density[i] < 0)
        	printf("water_vapour_density negative.\n.");
    }
    free(gravity_potential);
    for (int i = 0; i < NO_OF_LAYERS; ++i)
    {
        // horizontal wind fields are determind here
        for (int j = 0; j < NO_OF_VECTORS_H; ++j)
        {
            lat = latitude_vector[j];
            lon = longitude_vector[j];
            z_height = z_vector[NO_OF_SCALARS_H + j + i*NO_OF_VECTORS_PER_LAYER];
            // standard atmosphere: no wind
            if (TEST_ID == 0 || TEST_ID == 1)
                wind[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = 0;
            // JW test: specific wind field
            if (TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5 || TEST_ID == 6 || TEST_ID == 7)
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
                wind[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = u*cos(direction[j]);
            }
        }
    }
    for (int i = 0; i < NO_OF_LEVELS; ++i)
    {
    	#pragma omp parallel for private(lat, lon, z_height)
        for (int j = 0; j < NO_OF_SCALARS_H; ++j)
        {
            lat = latitude_scalar[j];
            lon = longitude_scalar[j];
            z_height = z_vector[j + i*NO_OF_VECTORS_PER_LAYER];
            if (TEST_ID == 0 || TEST_ID == 1 || TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5)
                wind[i*NO_OF_VECTORS_PER_LAYER + j] = 0;
        }
    }
    free(z_vector);
    free(latitude_scalar);
    free(longitude_scalar);
    free(direction);
    free(latitude_vector);
    free(longitude_vector);
    int scalar_dimid, vector_dimid, temp_id, density_dry_id, wind_id, density_vapour_id, density_liquid_id, density_solid_id, temperature_liquid_id, temperature_solid_id, ncid;
    if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "scalar_index", NO_OF_SCALARS, &scalar_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "vector_index", NO_OF_VECTORS, &vector_dimid)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperature_gas", NC_DOUBLE, 1, &scalar_dimid, &temp_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temp_id, "units", strlen("K"), "K")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "density_dry", NC_DOUBLE, 1, &scalar_dimid, &density_dry_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, density_dry_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "wind", NC_DOUBLE, 1, &vector_dimid, &wind_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, wind_id, "units", strlen("m/s"), "m/s")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "density_vapour", NC_DOUBLE, 1, &scalar_dimid, &density_vapour_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, density_vapour_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "density_liquid", NC_DOUBLE, 1, &scalar_dimid, &density_liquid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, density_liquid_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "density_solid", NC_DOUBLE, 1, &scalar_dimid, &density_solid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, density_solid_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperature_liquid", NC_DOUBLE, 1, &scalar_dimid, &temperature_liquid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temperature_liquid_id, "units", strlen("T"), "T")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperature_solid", NC_DOUBLE, 1, &scalar_dimid, &temperature_solid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temperature_solid_id, "units", strlen("T"), "T")))
        NCERR(retval);
    if ((retval = nc_enddef(ncid)))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temp_id, &temperature[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, density_dry_id, &rho[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, wind_id, &wind[0])))
        NCERR(retval);    
    if ((retval = nc_put_var_double(ncid, density_vapour_id, &water_vapour_density[0])))
        NCERR(retval);    
    if ((retval = nc_put_var_double(ncid, density_liquid_id, &liquid_water_density[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, density_solid_id, &solid_water_density[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temperature_liquid_id, &liquid_water_temp[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temperature_solid_id, &solid_water_temp[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
    	NCERR(retval);
    free(wind);
    free(pressure);
    free(temperature);
    free(rho);
    free(water_vapour_density);
    free(liquid_water_density);
    free(solid_water_density);
    free(liquid_water_temp);
    free(solid_water_temp);
    free(z_scalar);
    free(rel_humidity);
    free(OUTPUT_FILE);
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

















