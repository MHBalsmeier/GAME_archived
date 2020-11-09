/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
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
#include "../../shared/shared.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define P_0 100000.0
#define OMEGA (7.292115e-5)
#define C_D_P 1005.0
#define C_D_V 717.942189

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

// constants that are specific to the ICAO standard atmosphere
const double P_0_STANDARD = 101325;
const double TROPO_HEIGHT_STANDARD = 11e3;
const double INVERSE_HEIGHT_STANDARD = 20e3;
const double TEMP_GRADIENT_INV_STANDARD = 0.1/100;

int find_pressure_value(double, double, double *);
double sackur_tetrode(double, double);
double solve_specific_entropy_for_density(double, double);

int main(int argc, char *argv[])
{
	int TEST_ID;
   	TEST_ID = strtod(argv[1], NULL);
   	// determining the orography ID as a function of the test ID
	int ORO_ID;
	if (TEST_ID == 0 || TEST_ID == 8 || TEST_ID == 9)
	{
		ORO_ID = 0;
	}
	if (TEST_ID == 1)
	{
		ORO_ID = 1;
	}
	if (TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5)
	{
		ORO_ID = 2;
	}
	if (TEST_ID == 6 || TEST_ID == 7)
	{
		ORO_ID = 3;
	}
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
    const double TROPO_TEMP_STANDARD = T_SFC + TROPO_HEIGHT_STANDARD*TEMP_GRADIENT;
    double z_height;
    double lat, lon, u, v, eta, eta_v, T_perturb, distance, pressure_value, pressure_at_inv_standard, specific_humidity, total_density;
    double u_p = 1.0;
    double distance_scale = RADIUS/10;
    double lat_perturb = 2*M_PI/9;
    double lon_perturb = M_PI/9;
    // dummy arguments
    double dummy_0 = 0.0;
    double dummy_1 = 0.0;
    double dummy_2 = 0.0;
    double dummy_3 = 0.0;
    double dummy_4 = 0.0;
    double dummy_5 = 0.0;
    double dummy_6 = 0.0;
    int layer_index, h_index;
    int zero = 0;
    int one = 1;
    double one_double = 1;
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
            if (z_height < TROPO_HEIGHT_STANDARD)
            {
                temperature[i] = T_SFC + z_height*TEMP_GRADIENT;
                pressure[i] = P_0_STANDARD*pow(1 + TEMP_GRADIENT*z_height/T_SFC, -G/(R_D*TEMP_GRADIENT));
            }
            else if (z_height < INVERSE_HEIGHT_STANDARD)
            {
                temperature[i] = TROPO_TEMP_STANDARD;
                pressure[i] = P_0_STANDARD*pow(1 + TEMP_GRADIENT*TROPO_HEIGHT_STANDARD/T_SFC, -G/(R_D*TEMP_GRADIENT))*exp(-G*(z_height - TROPO_HEIGHT_STANDARD)/(R_D*TROPO_TEMP_STANDARD));
            }
            else
            {
            	temperature[i] = TROPO_TEMP_STANDARD + TEMP_GRADIENT_INV_STANDARD*(z_height - INVERSE_HEIGHT_STANDARD);
            	pressure_at_inv_standard = P_0_STANDARD*pow(1 + TEMP_GRADIENT*TROPO_HEIGHT_STANDARD/T_SFC, -G/(R_D*TEMP_GRADIENT))*exp(-G*(INVERSE_HEIGHT_STANDARD - TROPO_HEIGHT_STANDARD)/(R_D*TROPO_TEMP_STANDARD));
                pressure[i] = pressure_at_inv_standard*pow(1 + TEMP_GRADIENT*(z_height - INVERSE_HEIGHT_STANDARD)/T_SFC, -G/(R_D*TEMP_GRADIENT));
            }
        }
        // JW test
        if (TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5 || TEST_ID == 6 || TEST_ID == 7)
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
                {
                    rel_humidity[i] = 0.7;
                }
                else
                {
                    rel_humidity[i] = 0;
                }
            }
            else
            {
                temperature[i] = T_0*pow(eta, R_D*GAMMA/G) + DELTA_T*pow(ETA_T - eta, 5) + T_perturb;
            }
        }
        // dry Ullrich test
        if (TEST_ID == 8)
        {
        	baroclinic_wave_test(&one, &zero, &one, &one_double, &lon, &lat, &pressure[i], &z_height, &one, &dummy_0, &dummy_1, &temperature[i], &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
        }
        // moist Ullrich test
        if (TEST_ID == 9)
        {
        	baroclinic_wave_test(&one, &one, &one, &one_double, &lon, &lat, &pressure[i], &z_height, &one, &dummy_0, &dummy_1, &temperature[i], &dummy_2, &dummy_3, &dummy_4, &total_density, &specific_humidity);
        	water_vapour_density[i] = total_density*specific_humidity;
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
        	lower_entropy_value = sackur_tetrode(rho[i + NO_OF_SCALARS_H], temperature[i + NO_OF_SCALARS_H]);
        	temperature_mean = 0.5*(temperature[i] + temperature[i + NO_OF_SCALARS_H]);
        	delta_temperature = temperature[i] - temperature[i + NO_OF_SCALARS_H];
        	delta_gravity_potential = gravity_potential[i] - gravity_potential[i + NO_OF_SCALARS_H];
        	entropy_value = lower_entropy_value + (delta_gravity_potential + C_D_P*delta_temperature)/temperature_mean;
        	rho[i] = solve_specific_entropy_for_density(entropy_value, temperature[i]);
        }
        if (TEST_ID == 4 || TEST_ID == 5 || TEST_ID == 7)
        {
		    water_vapour_density[i] = water_vapour_density_from_rel_humidity(rel_humidity[i], temperature[i], rho[i]);
		    if (water_vapour_density[i] < 0)
		    {
		    	printf("water_vapour_density negative.\n.");
			}
    	}
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
            {
                wind[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = 0;
            }
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
            // dry Ullrich test
            if (TEST_ID == 8)
            {
        		baroclinic_wave_test(&one, &zero, &one, &one_double, &lon, &lat, &dummy_0, &z_height, &one, &u, &v, &dummy_1, &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
                wind[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = u*cos(direction[j]) + v*sin(direction[j]);
            }
            // moist Ullrich test
            if (TEST_ID == 9)
            {
        		baroclinic_wave_test(&one, &one, &one, &one_double, &lon, &lat, &dummy_0, &z_height, &one, &u, &v, &dummy_1, &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
                wind[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = u*cos(direction[j]) + v*sin(direction[j]);
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
            if (TEST_ID == 0 || TEST_ID == 1 || TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5 || TEST_ID == 6 || TEST_ID == 7 || TEST_ID == 8 || TEST_ID == 9)
            {
                wind[i*NO_OF_VECTORS_PER_LAYER + j] = 0;
            }
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

double sackur_tetrode(double mass_density, double temperature)
{
	double mean_particle_mass = 0.004810e-23;
	double entropy_constant = 2429487178047751925300627872548148580712448.000000;
	double particle_density = mass_density/mean_particle_mass;
	// returns the specific entropy as a function of the mass density and the temperature
	double result;
    result = K_B*(3.0/2*log(entropy_constant)
    + log(1/particle_density)
    + 3.0/2*log(mean_particle_mass*C_D_V*temperature))/mean_particle_mass;
    return result;
}

double solve_specific_entropy_for_density(double specific_entropy, double temperature)
{
	// returns the density as a function of the specific entropy and the temperature
	double mean_particle_mass = 0.004810e-23;
	double entropy_constant = 2429487178047751925300627872548148580712448.000000;
    double particle_density = exp(-mean_particle_mass*specific_entropy/K_B)*pow(entropy_constant, 1.5)*pow(mean_particle_mass*C_D_V*temperature, 1.5);
    double result = particle_density*mean_particle_mass;
    return result;
}











