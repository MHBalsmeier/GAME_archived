/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, the initial state of the simulation is read in from a netcdf file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include <atmostracers.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../spatial_operators/spatial_operators.h"
#include "../thermodynamics/thermodynamics.h"
#include "../../grid_generator/src/standard.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

int set_soil_temp(Grid *grid, Soil *, State *, double []);

int set_ideal_init(State *state, Grid* grid, Dualgrid* dualgrid, Soil *soil, int test_id, char geo_prop_file[])
{
	/*
	This function sets the initial state of the model atmosphere for idealized test cases.
	*/
	
    // reading the grid properties which are not part of the struct grid
    double *latitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *longitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    int ncid_grid, retval, latitude_vector_id, longitude_vector_id;
    if ((retval = nc_open(geo_prop_file, NC_NOWRITE, &ncid_grid)))
        NCERR(retval);
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
	
    double *pressure = malloc(NO_OF_SCALARS*sizeof(double));
    double *temperature = malloc(NO_OF_SCALARS*sizeof(double));
    double *water_vapour_density = malloc(NO_OF_SCALARS*sizeof(double));
    double z_height;
    double lat, lon, u, v, pressure_value, specific_humidity, total_density;
    // dummy argument
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
    #pragma omp parallel for private(layer_index, h_index, lat, lon, z_height)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
        lat = grid -> latitude_scalar[h_index];
        lon = grid -> longitude_scalar[h_index];
        z_height = grid -> z_scalar[i];
        // standard atmosphere
        if (test_id == 0 || test_id == 1 || test_id == 2)
        {
            temperature[i] = standard_temp(z_height);
            pressure[i] = standard_pres(z_height);
        }
        // dry Ullrich test
        if (test_id == 3 || test_id == 4 || test_id == 5)
        {
        	baroclinic_wave_test(&one, &zero, &one, &one_double, &lon, &lat, &pressure[i], &z_height, &one, &dummy_0, &dummy_1, &temperature[i],
        	&dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
        }
        // moist Ullrich test
        if (test_id == 6 || test_id == 7 || test_id == 8)
        {
        	baroclinic_wave_test(&one, &one, &one, &one_double, &lon, &lat, &pressure[i], &z_height, &one, &dummy_0, &dummy_1, &temperature[i],
        	&dummy_2, &dummy_4, &dummy_5, &total_density, &specific_humidity);
        	water_vapour_density[i] = total_density*specific_humidity;
        }
    }

    // horizontal wind fields are determind here
    #pragma omp parallel for private(lat, lon, z_height, u, v, dummy_0, dummy_1, dummy_2, dummy_3, dummy_4, dummy_5, dummy_6)
    for (int i = 0; i < NO_OF_LAYERS; ++i)
    {
        for (int j = 0; j < NO_OF_VECTORS_H; ++j)
        {
            lat = latitude_vector[j];
            lon = longitude_vector[j];
            z_height = grid -> z_vector[NO_OF_SCALARS_H + j + i*NO_OF_VECTORS_PER_LAYER];
            // standard atmosphere: no wind
            if (test_id == 0 || test_id == 1 || test_id == 2)
            {
                state -> wind[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = 0;
            }
            // dry Ullrich test
            if (test_id == 3 || test_id == 4 || test_id == 5)
            {
        		baroclinic_wave_test(&one, &zero, &one, &one_double, &lon, &lat, &dummy_0, &z_height, &one, &u, &v, &dummy_1, &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
                state -> wind[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = u*cos(grid -> direction[j]) + v*sin(grid -> direction[j]);
            }
            // moist Ullrich test
            if (test_id == 6 || test_id == 7 || test_id == 8)
            {
        		baroclinic_wave_test(&one, &one, &one, &one_double, &lon, &lat, &dummy_0, &z_height, &one, &u, &v, &dummy_1, &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
                state -> wind[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = u*cos(grid -> direction[j]) + v*sin(grid -> direction[j]);
            }
        }
    }
    free(latitude_vector);
    free(longitude_vector);
    // setting the vertical wind field equal to zero
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_LEVELS; ++i)
    {
        for (int j = 0; j < NO_OF_SCALARS_H; ++j)
        {
            state -> wind[i*NO_OF_VECTORS_PER_LAYER + j] = 0;
        }
    }    
    
    Diagnostics *diagnostics = calloc(1, sizeof(Diagnostics));
    Forcings *forcings = calloc(1, sizeof(Forcings));
    // this is the density which has not yet been hydrostatically balanced
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = pressure[i]/(specific_gas_constants(0)*temperature[i]);
	}
	scalar_times_vector(diagnostics -> scalar_field_placeholder, state -> wind, diagnostics -> flux_density, grid);
	// Now, the potential vorticity is evaluated.
	calc_pot_vort(state -> wind, diagnostics -> scalar_field_placeholder, diagnostics, grid, dualgrid);
	// Now, the generalized Coriolis term is evaluated.
	vorticity_flux(diagnostics -> flux_density, diagnostics -> pot_vort, forcings -> pot_vort_tend, grid, dualgrid);
	
	// Kinetic energy is prepared for the gradient term of the Lamb transformation.
	inner_product(state -> wind, state -> wind, diagnostics -> v_squared, grid);
    // density is determined out of the hydrostatic equation
    int scalar_index;
    double b, c;
    // theta_pert and exner_pert are a misuse of name here, they contain the full values here
    #pragma omp parallel for private(scalar_index, b, c, pressure_value)
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		// integrating from bottom to top
		for (int layer_index = NO_OF_LAYERS - 1; layer_index >= 0; --layer_index)
		{
			scalar_index = layer_index*NO_OF_SCALARS_H + h_index;
			// lowest layer
			if (layer_index == NO_OF_LAYERS - 1)
			{
				pressure_value = pressure[scalar_index];
				state -> exner_pert[scalar_index] = pow(pressure_value/P_0, specific_gas_constants(0)/spec_heat_capacities_p_gas(0));
				state -> theta_pert[scalar_index] = temperature[scalar_index]/state -> exner_pert[scalar_index];
			}
			// other layers
			else
			{
				// solving a quadratic equation for the Exner pressure
				b = -0.5*state -> exner_pert[scalar_index + NO_OF_SCALARS_H]/temperature[scalar_index + NO_OF_SCALARS_H]
				*(temperature[scalar_index] - temperature[scalar_index + NO_OF_SCALARS_H]
				+ 2/spec_heat_capacities_p_gas(0)*(grid -> gravity_potential[scalar_index] - grid -> gravity_potential[scalar_index + NO_OF_SCALARS_H]
				+ 0.5*diagnostics -> v_squared[scalar_index] - 0.5*diagnostics -> v_squared[scalar_index + NO_OF_SCALARS_H]
				- (grid -> z_scalar[scalar_index] - grid -> z_scalar[scalar_index + NO_OF_SCALARS_H])*forcings -> pot_vort_tend[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER]));
				c = pow(state -> exner_pert[scalar_index + NO_OF_SCALARS_H], 2)*temperature[scalar_index]/temperature[scalar_index + NO_OF_SCALARS_H];
				state -> exner_pert[scalar_index] = b + pow((pow(b, 2) + c), 0.5);
			}
			// scalar_field_placeholder is the dry air density here
			diagnostics -> scalar_field_placeholder[scalar_index] = P_0*pow(state -> exner_pert[scalar_index],
			spec_heat_capacities_p_gas(0)/specific_gas_constants(0))/(specific_gas_constants(0)*temperature[scalar_index]);
		}
	}
    free(pressure);
    free(forcings);
    
    double *temperatures = malloc((NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS*sizeof(double));
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		for (int j = 0; j < NO_OF_CONDENSED_CONSTITUENTS; ++j)
		{
			// condensed densities are zero in all test states
			state -> rho[j*NO_OF_SCALARS + i] = 0;
			// a local LTE is assumed in all test states
			temperatures[j*NO_OF_SCALARS + i] = temperature[i];
		}
		// the dry air density
		state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i] = diagnostics -> scalar_field_placeholder[i];
		// water vapout density
		if (NO_OF_CONDENSED_CONSTITUENTS == 4)
		{
			state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i] = water_vapour_density[i];
		}
		// gas temperature
		temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i] = temperature[i];
	}
	free(diagnostics);
    free(temperature);
    free(water_vapour_density);
    
    // setting the soil temperature
    set_soil_temp(grid, soil, state, temperatures);
}

int read_init_data(char FILE_NAME[], State *state, Grid* grid, Soil *soil)
{
	/*
	This function reads the initial state of the model atmosphere from a netcdf file.
	*/
	
    double *temperatures = malloc((NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS*sizeof(double));
    double *sst = malloc(NO_OF_SCALARS_H*sizeof(double));
    int retval, ncid;
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        NCERR(retval);
    int densities_id, temperatures_id, wind_id, sst_id, stretching_parameter_id, sst_avail;
    double stretching_parameter;
    if ((retval = nc_inq_varid(ncid, "densities", &densities_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperatures", &temperatures_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "wind", &wind_id)))
        NCERR(retval);
    // figuring out if the netcdf file contains SST
    sst_avail = 0;
    if (nc_inq_varid(ncid, "sst", &sst_id) == 0)
    {
    	sst_avail = 1;
    	printf("SST found in initialization file.\n");
    }
    else
    {	
    	printf("SST not found in initialization file.\n");
    }
    if ((retval = nc_inq_varid(ncid, "stretching_parameter", &stretching_parameter_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, densities_id, &state -> rho[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, temperatures_id, &temperatures[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, wind_id, &state -> wind[0])))
        NCERR(retval);
    // reading the SST data if it is present in the netcdf file
    if (sst_avail == 1)
    {
		if ((retval = nc_get_var_double(ncid, sst_id, &sst[0])))
		    NCERR(retval);
    }
    if ((retval = nc_get_var_double(ncid, stretching_parameter_id, &stretching_parameter)))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
        NCERR(retval);
    
	// resricting the maximum relative humidity to 100 %
    if (NO_OF_CONDENSED_CONSTITUENTS == 4)
    {
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			if (rel_humidity(state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i], temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]) > 1)
			{
				state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i] = state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i]
				/rel_humidity(state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i], temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]);
			}
		}
    }
	
    // determining the temperature densities of the condensates
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS; ++i)
	{
		state -> condensed_density_temperatures[i] = state -> rho[i]*temperatures[i];
	}
	
	// diagnostic thermodynamical quantities
	double pressure, pot_temp;
	#pragma omp parallel for private(pressure, pot_temp)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		pressure = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*specific_gas_constants(0)*temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i];
		pot_temp = temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*pow(P_0/pressure, specific_gas_constants(0)/spec_heat_capacities_p_gas(0));
		state -> rhotheta[i] = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*pot_temp;
		// calculating the potential temperature perturbation
		state -> theta_pert[i] = pot_temp - grid -> theta_bg[i];
		// calculating the Exner pressure perturbation
		state -> exner_pert[i] = temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]/(grid -> theta_bg[i] + state -> theta_pert[i]) - grid -> exner_bg[i];
	}
	
    // checks
    // checking for negative densities
    # pragma omp parallel for
	for (int i = 0; i < NO_OF_CONSTITUENTS*NO_OF_SCALARS; ++i)
	{
		if (state -> rho[i] < 0)
	    {
			printf("Negative density found.\n");
			printf("Aborting.\n");
			exit(1);
    	}
	}
	// checking wether the stretching parameters of the grid used for creating the input file and the grid file read in conform
    if (grid -> stretching_parameter != stretching_parameter)
    {
    	printf("Stretching parameters of grid and input file do not conform.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
    
    free(sst);
    free(temperatures);
    
    // setting the soil temperature
    set_soil_temp(grid, soil, state, temperatures);
    
    return 0;
}

int set_soil_temp(Grid *grid, Soil *soil, State *state, double temperatures[])
{
	/*
	This function sets the soil temperature.
	*/
	
	int soil_index;
	double z_soil, t_sfc;
	// setting the soil temperature
	#pragma omp parallel for private(soil_index, z_soil, t_sfc)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// temperature at the surface
		// land surface or sea surface if SST is not available
		if (grid -> is_land[i] == 1)
		{
			t_sfc = temperatures[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + NO_OF_SCALARS - NO_OF_SCALARS_H + i];
		}
		// sea surface if SST is available
		else
		{
			t_sfc = grid -> t_const_soil;
		}
		
		// loop over all soil layers
		for (int soil_layer_index = 0; soil_layer_index < NO_OF_SOIL_LAYERS; ++soil_layer_index)
		{
			// index of this soil grid point
			soil_index = i + soil_layer_index*NO_OF_SCALARS_H;
			z_soil = grid -> z_t_const/NO_OF_SOIL_LAYERS*(0.5 + soil_layer_index);
			soil -> temperature[soil_index] = t_sfc + (grid -> t_const_soil - t_sfc)*z_soil/grid -> z_t_const;
		}
	}
	// returning 0 indicating success
	return 0;
}



