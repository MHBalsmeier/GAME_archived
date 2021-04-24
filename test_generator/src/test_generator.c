/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include <stdlib.h>
#include "../../core/src/enum_and_typedefs.h"
#include "../../core/src/thermodynamics/thermodynamics.h"
#include "../../core/src/settings.h"
#include "../../core/src/spatial_operators/spatial_operators.h"
#include "../../core/src/io/io.h"
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

int main(int argc, char *argv[])
{
	// some thermodynamical quantities
	double R_D = specific_gas_constants_lookup(0);
	double C_D_P = spec_heat_capacities_p_gas_lookup(0);
	int TEST_ID;
   	TEST_ID = strtod(argv[1], NULL);
   	int NO_OF_ORO_LAYERS = strtod(argv[2], NULL);
   	const int VERT_GRID_TYPE = strtod(argv[3], NULL);
   	const double TOA = strtof(argv[4], NULL);
   	if (VERT_GRID_TYPE == 1)
   	{
   		NO_OF_ORO_LAYERS = 0;
   	}
   	
   	// determining the orography ID as a function of the test ID
	int ORO_ID;
	if (TEST_ID == 0 || TEST_ID == 8 || TEST_ID == 9)
	{
		ORO_ID = 0;
	}
	if (TEST_ID == 1 || TEST_ID == 13)
	{
		ORO_ID = 1;
	}
	if (TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5)
	{
		ORO_ID = 2;
	}
	if (TEST_ID == 6 || TEST_ID == 7 || TEST_ID == 10 || TEST_ID == 11 || TEST_ID == 12 || TEST_ID == 14)
	{
		ORO_ID = 3;
	}
    int FILE_NAME_LENGTH = 100;
    char *GEO_PROP_FILE_PRE = malloc((FILE_NAME_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE_PRE, "../grid_generator/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    FILE_NAME_LENGTH = strlen(GEO_PROP_FILE_PRE);
    free(GEO_PROP_FILE_PRE);
    char *GEO_PROP_FILE = malloc((FILE_NAME_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE, "../grid_generator/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    Grid *grid = calloc(1, sizeof(Grid));
    Dualgrid *dualgrid = calloc(1, sizeof(Dualgrid));
    set_grid_properties(grid, dualgrid, GEO_PROP_FILE);
    
    int OUTPUT_FILE_LENGTH = 100;
    char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE_PRE, "test_states/test_%d_B%dL%dT%d_O%d_OL%d_SCVT.nc", TEST_ID, RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((FILE_NAME_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "test_states/test_%d_B%dL%dT%d_O%d_OL%d_SCVT.nc", TEST_ID, RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    double *pressure = malloc(NO_OF_SCALARS*sizeof(double));
    double *temperature = malloc(NO_OF_SCALARS*sizeof(double));
    double *rel_humidity = malloc(NO_OF_SCALARS*sizeof(double));
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
    State *state = calloc(1, sizeof(State));
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
        lat = grid -> latitude_scalar[h_index];
        lon = grid -> longitude_scalar[h_index];
        z_height = grid -> z_scalar[i];
        rel_humidity[i] = 0;
        // standard atmosphere
        if (TEST_ID == 0 || TEST_ID == 1 || TEST_ID == 12)
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
        if (TEST_ID == 8 || TEST_ID == 10 || TEST_ID == 13 || TEST_ID == 14)
        {
        	baroclinic_wave_test(&one, &zero, &one, &one_double, &lon, &lat, &pressure[i], &z_height, &one, &dummy_0, &dummy_1, &temperature[i], &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
        }
        // moist Ullrich test
        if (TEST_ID == 9 || TEST_ID == 11)
        {
        	baroclinic_wave_test(&one, &one, &one, &one_double, &lon, &lat, &pressure[i], &z_height, &one, &dummy_0, &dummy_1, &temperature[i], &dummy_2, &dummy_3, &dummy_4, &total_density, &specific_humidity);
        	water_vapour_density[i] = total_density*specific_humidity;
        }
	    liquid_water_density[i] = 0;
	    solid_water_density[i] = 0;
        liquid_water_temp[i] = temperature[i];
        solid_water_temp[i] = temperature[i];
    }
    Forcings *forcings = calloc(1, sizeof(Forcings));
    
    // reading the grid properties which are not part of the struct grid
    double *latitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *longitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    int ncid_grid, retval, latitude_vector_id, longitude_vector_id;
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid_grid)))
        NCERR(retval);
    free(GEO_PROP_FILE);
    if ((retval = nc_inq_varid(ncid_grid, "latitude_vector", &latitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_vector", &longitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitude_vector_id, &latitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitude_vector_id, &longitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid_grid)))
        NCERR(retval);

    // horizontal wind fields are determind here
    for (int i = 0; i < NO_OF_LAYERS; ++i)
    {
        for (int j = 0; j < NO_OF_VECTORS_H; ++j)
        {
            lat = latitude_vector[j];
            lon = longitude_vector[j];
            z_height = grid -> z_vector[NO_OF_SCALARS_H + j + i*NO_OF_VECTORS_PER_LAYER];
            // standard atmosphere: no wind
            if (TEST_ID == 0 || TEST_ID == 1 || TEST_ID == 12)
            {
                state -> velocity_gas[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = 0;
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
                state -> velocity_gas[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = u*cos(grid -> direction[j]);
            }
            // dry Ullrich test
            if (TEST_ID == 8 || TEST_ID == 10 || TEST_ID == 13)
            {
        		baroclinic_wave_test(&one, &zero, &one, &one_double, &lon, &lat, &dummy_0, &z_height, &one, &u, &v, &dummy_1, &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
                state -> velocity_gas[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = u*cos(grid -> direction[j]) + v*sin(grid -> direction[j]);
            }
            if (TEST_ID == 14)
            {
        		baroclinic_wave_test(&one, &zero, &one, &one_double, &lon, &lat, &dummy_0, &z_height, &one, &u, &v, &dummy_1, &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
                state -> velocity_gas[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = -(u*cos(grid -> direction[j]) + v*sin(grid -> direction[j]));
            }
            // moist Ullrich test
            if (TEST_ID == 9 || TEST_ID == 11)
            {
        		baroclinic_wave_test(&one, &one, &one, &one_double, &lon, &lat, &dummy_0, &z_height, &one, &u, &v, &dummy_1, &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
                state -> velocity_gas[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = u*cos(grid -> direction[j]) + v*sin(grid -> direction[j]);
            }
        }
    }
    // setting the vertical wind field equal to zero
    for (int i = 0; i < NO_OF_LEVELS; ++i)
    {
    	#pragma omp parallel for private(lat, lon, z_height)
        for (int j = 0; j < NO_OF_SCALARS_H; ++j)
        {
            lat = grid -> latitude_scalar[j];
            lon = grid -> longitude_scalar[j];
            z_height = grid -> z_vector[j + i*NO_OF_VECTORS_PER_LAYER];
            if (TEST_ID == 0 || TEST_ID == 1 || TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5 || TEST_ID == 6 || TEST_ID == 7 || TEST_ID == 8 || TEST_ID == 9 || TEST_ID == 10 || TEST_ID == 11 || TEST_ID == 12 || TEST_ID == 13 || TEST_ID == 14)
            {
                state -> velocity_gas[i*NO_OF_VECTORS_PER_LAYER + j] = 0;
            }
        }
    }    
    
    Diagnostics *diagnostics = calloc(1, sizeof(Diagnostics));
    // this is the density which has not yet been hydrostatically balanced
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = pressure[i]/(R_D*temperature[i]);
	}
	scalar_times_vector(diagnostics -> scalar_field_placeholder, state -> velocity_gas, diagnostics -> flux_density, grid);
	// Now, the potential vorticity is evaluated.
	calc_pot_vort(state -> velocity_gas, diagnostics -> scalar_field_placeholder, diagnostics, grid, dualgrid);
	// Now, the generalized Coriolis term is evaluated.
	vorticity_flux(diagnostics -> flux_density, diagnostics -> pot_vort, forcings -> pot_vort_tend, grid, dualgrid);
	free(dualgrid);
	// Kinetic energy is prepared for the gradient term of the Lamb transformation.
	inner_product(state -> velocity_gas, state -> velocity_gas, diagnostics -> e_kin, grid);
	// Taking the gradient of the kinetic energy
	grad(diagnostics -> e_kin, forcings -> e_kin_grad, grid);
    // density is determined out of the hydrostatic equation
    double entropy_value, temperature_mean, delta_temperature, delta_gravity_potential, lower_entropy_value, delta_z;
    for (int i = NO_OF_SCALARS - 1; i >= 0; --i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
    	// at the lowest layer the density is set using the equation of state, can be considered a boundary condition
    	if (layer_index == NO_OF_LAYERS - 1)
    	{
        	state -> mass_densities[i] = pressure[i]/(R_D*temperature[i]);
        }
        else
        {
        	lower_entropy_value = spec_entropy_from_temp(state -> mass_densities[i + NO_OF_SCALARS_H], temperature[i + NO_OF_SCALARS_H]);
        	temperature_mean = 0.5*(temperature[i] + temperature[i + NO_OF_SCALARS_H]);
        	delta_temperature = temperature[i] - temperature[i + NO_OF_SCALARS_H];
        	delta_z = grid -> z_scalar[i] - grid -> z_scalar[i + NO_OF_SCALARS_H];
        	delta_gravity_potential = grid -> gravity_potential[i] - grid -> gravity_potential[i + NO_OF_SCALARS_H]
        	- delta_z*(forcings -> pot_vort_tend[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + h_index] - 0.5*forcings -> e_kin_grad[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + h_index]);
        	entropy_value = lower_entropy_value + (delta_gravity_potential + C_D_P*delta_temperature)/temperature_mean;
        	state -> mass_densities[i] = solve_specific_entropy_for_density(entropy_value, temperature[i]);
        }
        if (TEST_ID == 4 || TEST_ID == 5 || TEST_ID == 7)
        {
		    water_vapour_density[i] = water_vapour_density_from_rel_humidity(rel_humidity[i], temperature[i], state -> mass_densities[i]);
		    if (water_vapour_density[i] < 0)
		    {
		    	printf("water_vapour_density negative.\n.");
			}
    	}
    }
    
    free(latitude_vector);
    free(longitude_vector);
    free(state);
    free(forcings);
    int scalar_dimid, vector_dimid, temp_id, density_dry_id, wind_id, density_vapour_id, density_liquid_id, density_solid_id, temperature_liquid_id, temperature_solid_id, ncid, single_double_dimid, stretching_parameter_id;
    if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "scalar_index", NO_OF_SCALARS, &scalar_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "vector_index", NO_OF_VECTORS, &vector_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "single_double_dimid_index", 1, &single_double_dimid)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "stretching_parameter", NC_DOUBLE, 1, &single_double_dimid, &stretching_parameter_id)))
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
    if ((retval = nc_put_var_double(ncid, stretching_parameter_id, &grid -> stretching_parameter)))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temp_id, &temperature[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, density_dry_id, &state -> mass_densities[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, wind_id, &state -> velocity_gas[0])))
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
    free(grid);
    free(pressure);
    free(temperature);
    free(water_vapour_density);
    free(liquid_water_density);
    free(solid_water_density);
    free(liquid_water_temp);
    free(solid_water_temp);
    free(rel_humidity);
    free(OUTPUT_FILE);
    return 0;
}

int find_pressure_value(double lat, double z_height, double *result)
{
	// this function finds the pressure at a given height (as a function of latitude) for the JW test by iterative calls to the function find_z_from_p_jw
    double p = P_0/2;
    double precision = 0.0001;
    double z;
    double current_max = P_0;
    double current_min = 0;
    find_z_from_p_jw(lat, p, &z);
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
        find_z_from_p_jw(lat, p, &z);
    }
    *result = p;
    return 0;
}






