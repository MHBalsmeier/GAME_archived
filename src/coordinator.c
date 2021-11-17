/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
The main organizes the model, manages the time stepping, calls model output, collects the lowest model layer wind for 10 m wind mean and so on. All the memory needed for the integration is allocated and freed here.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <geos95.h>
#include "game_types.h"
#include "game_constants.h"
#include "io/io.h"
#include "spatial_operators/spatial_operators.h"
#include "radiation/radiation.h"
#include "thermodynamics.h"
#include "time_stepping/time_stepping.h"
#include "soil.h"

int main(int argc, char *argv[])
{
    // taking the timestamp to measure the performance
    clock_t begin = clock();
    
    /*
    Reading command line input.
    ---------------------------
    */
    Grid *grid = calloc(1, sizeof(Grid));
    Dualgrid *dualgrid = calloc(1, sizeof(Dualgrid));
    size_t len = strlen(argv[1]);
    char *TOTAL_RUN_SPAN_PRE = malloc((len + 1)*sizeof(char));
    strcpy(TOTAL_RUN_SPAN_PRE, argv[1]);
    int TOTAL_RUN_SPAN = strtol(TOTAL_RUN_SPAN_PRE, NULL, 10);
    free(TOTAL_RUN_SPAN_PRE);
    len = strlen(argv[2]);
    char *WRITE_OUT_INTERVAL_PRE = malloc((len + 1)*sizeof(char));
    strcpy(WRITE_OUT_INTERVAL_PRE, argv[2]);
    int WRITE_OUT_INTERVAL = strtol(WRITE_OUT_INTERVAL_PRE, NULL, 10);
    free(WRITE_OUT_INTERVAL_PRE);
    double dt_parameter = strtod(argv[3], NULL);
    Config *config = calloc(1, sizeof(Config));
    config -> momentum_diff_h = strtod(argv[4], NULL);
    config -> momentum_diff_v = strtod(argv[5], NULL);
    config -> rad_on = strtod(argv[6], NULL);
    int write_out_mass_integrals;
    write_out_mass_integrals = strtod(argv[7], NULL);
    int write_out_rhotheta_integral; 
    write_out_rhotheta_integral = strtod(argv[8], NULL);
    int write_out_energy_integrals;
    write_out_energy_integrals = strtod(argv[9], NULL);
    config -> temperature_diff_h = strtod(argv[10], NULL);
    double radiation_delta_t;
    radiation_delta_t = strtod(argv[11], NULL);
    int year;
    year = strtod(argv[12], NULL);
    int month;
    month = strtod(argv[13], NULL);
    len = strlen(argv[13]);
    char *month_string = malloc((len + 1)*sizeof(char));
    strcpy(month_string, argv[13]);
    int day;
    day = strtod(argv[14], NULL);
    len = strlen(argv[14]);
    char *day_string = malloc((len + 1)*sizeof(char));
    strcpy(day_string, argv[14]);
    int hour;
    hour = strtod(argv[15], NULL);
    len = strlen(argv[15]);
    char *hour_string = malloc((len + 1)*sizeof(char));
    strcpy(hour_string, argv[15]);
    config -> temperature_diff_v = strtod(argv[16], NULL);
    len = strlen(argv[17]);
    char *RUN_ID = malloc((len + 1)*sizeof(char));
    strcpy(RUN_ID, argv[17]);
    int toa;
	toa = strtod(argv[18], NULL);
	int ORO_ID;
	ORO_ID = strtod(argv[19], NULL);
    int IDEAL_INPUT_ID;
    IDEAL_INPUT_ID = strtod(argv[20], NULL);
    Config_io *config_io = calloc(1, sizeof(Config_io));
	config_io -> grib_output_switch = strtod(argv[21], NULL);
	config_io -> netcdf_output_switch = strtod(argv[22], NULL);
	config_io -> pressure_level_output_switch = strtod(argv[23], NULL);
	config_io -> model_level_output_switch = strtod(argv[24], NULL);
	config_io -> surface_output_switch = strtod(argv[25], NULL);
	grid -> no_of_oro_layers = strtod(argv[26], NULL);
	int VERT_GRID_TYPE = strtod(argv[27], NULL);
	config -> assume_lte = strtod(argv[28], NULL);
	config -> slow_fast_ratio = strtod(argv[29], NULL);
	config -> delta_t_between_analyses = strtod(argv[30], NULL);
	config -> diff_h_smag_fac = strtod(argv[31], NULL);
	config -> shear_bg = strtod(argv[32], NULL);
	config -> damping_start_height_over_toa = strtod(argv[33], NULL);
	config -> damping_coeff_max = strtod(argv[34], NULL);
	config -> explicit_boundary_layer = strtod(argv[35], NULL);
	config -> tracer_diff_h = strtod(argv[36], NULL);
	config -> tracer_diff_v = strtod(argv[37], NULL);
	config -> impl_thermo_weight = strtod(argv[38], NULL);
	config -> cloud_droplets_velocity = strtod(argv[39], NULL);
	config -> precipitation_droplets_velocity = strtod(argv[40], NULL);
	
	/*
	Checking user input for correctness:
	------------------------------------
	*/
    if (grid -> no_of_oro_layers < 0 || grid -> no_of_oro_layers >= NO_OF_LAYERS)
    {
    	printf("It must be 0 <= orography_layers < NO_OF_LAYERS.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
    if (WRITE_OUT_INTERVAL < 900)
    {
    	printf("It is WRITE_OUT_INTERVAL < 900.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
	if (config -> assume_lte != 0 && config -> assume_lte != 1)
	{
		printf("simplified_moisture_switch must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	// in the case of block-shaped mountains, no layers follow the orography
	if (VERT_GRID_TYPE == 1)
	{
		grid -> no_of_oro_layers = 0;
	}
	if (config_io -> grib_output_switch == 0 && config_io -> netcdf_output_switch == 0)
	{
		printf("Either grib_output_switch or netcdf_output_switch must be set to 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> momentum_diff_h == 0 && config -> momentum_diff_v == 1)
	{
		printf("Horizontal momentum diffusion cannot be off if vertical momentum diffusion is on.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> temperature_diff_h == 0 && config -> temperature_diff_v == 1)
	{
		printf("Horizontal temperature diffusion cannot be off if vertical temperature diffusion is on.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> tracer_diff_h == 0 && config -> tracer_diff_v == 1)
	{
		printf("Horizontal tracer diffusion cannot be off if vertical tracer diffusion is on.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	
    /*
    This sets the ORO_ID (orography ID) as a function of the IDEAL_INPUT_ID.
    ---------------------------------------------------------------------------
    */
	if (IDEAL_INPUT_ID == 0 || IDEAL_INPUT_ID == 3 || IDEAL_INPUT_ID == 6)
    {
		ORO_ID = 0;
    }
	if (IDEAL_INPUT_ID == 1 || IDEAL_INPUT_ID == 4 || IDEAL_INPUT_ID == 7)
    {
		ORO_ID = 1;
    }
	if (IDEAL_INPUT_ID == 2 || IDEAL_INPUT_ID == 5 || IDEAL_INPUT_ID == 8)
    {
		ORO_ID = 2;
    }
    /*
	Determining the name of the grid file from the RES_ID, NO_OF_LAYERS and so on.
    ------------------------------------------------------------------------------
    */
    char GEO_PROP_FILE_PRE[200];
	sprintf(GEO_PROP_FILE_PRE, "../../grid_generator/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NO_OF_LAYERS, toa, ORO_ID, grid -> no_of_oro_layers);
    char GEO_PROP_FILE[strlen(GEO_PROP_FILE_PRE) + 1];
    strcpy(GEO_PROP_FILE, GEO_PROP_FILE_PRE);
    char SFC_PROP_FILE_PRE[200];
	sprintf(SFC_PROP_FILE_PRE, "../../surface_generator/surface_files/B%d_O%d_SCVT.nc", RES_ID, 2);
    char SFC_PROP_FILE[strlen(SFC_PROP_FILE_PRE) + 1];
    strcpy(SFC_PROP_FILE, SFC_PROP_FILE_PRE);
    
	// Determining the name of the init state file from the IDEAL_INPUT_ID, RES_ID, NO_OF_LAYERS and so on.
    char INIT_STATE_FILE_PRE[200];
    // The NWP case.
    if (IDEAL_INPUT_ID == -1)
    {
    	config -> nwp_mode = 1;
    	sprintf(INIT_STATE_FILE_PRE, "../../nwp_init/%d%s%s%s_B%dL%dT%d_O%d_OL%d_SCVT.nc", year, month_string, day_string, hour_string, RES_ID, NO_OF_LAYERS, toa, ORO_ID, grid -> no_of_oro_layers);
    }
    // The idealized input case.
    else
    {
    	config -> nwp_mode = 0;
		sprintf(INIT_STATE_FILE_PRE, "../../test_generator/test_states/test_%d_B%dL%dT%d_O%d_OL%d_SCVT.nc", IDEAL_INPUT_ID, RES_ID, NO_OF_LAYERS, toa, ORO_ID, grid -> no_of_oro_layers);
    }
    char INIT_STATE_FILE[strlen(INIT_STATE_FILE_PRE) + 1];
    strcpy(INIT_STATE_FILE, INIT_STATE_FILE_PRE);
    
    /*
    Determining the Unix time stamp of the initialization (UTC).
    ------------------------------------------------------------
    */
    struct tm init_tm;
    init_tm.tm_year = year - 1900;
    init_tm.tm_mon = month - 1;
    init_tm.tm_mday = day;
    init_tm.tm_hour = hour;
    init_tm.tm_min = 0;
    init_tm.tm_sec = 0;
    // turning off DST
    init_tm.tm_isdst = 0;
    time_t init_time = mktime(&init_tm);
    // converting to double in UTC
    double t_init = (double) init_time + init_tm.tm_gmtoff;
    
    /*
    Giving the user some information on the run to about to be executed.
    --------------------------------------------------------------------
    */
    char *stars  = malloc(83*sizeof(char));
    for (int i = 0; i < 81; ++i)
        stars[i] = '*';
    stars[81] = '\n';
    stars[82] = '\0';
    printf("%s", stars);
    printf("*\t\t\t\t\t\t\t\t\t\t*\n");
    printf("*\t\t\t\tThis is the GAME\t\t\t\t*\n");
    printf("*\t\t\tGeophysical Fluids Modeling Framework\t\t\t*\n");
    printf("*\t\t\t\t\t\t\t\t\t\t*\n");
    printf("%s", stars);
    printf("Released under the MIT license, visit https://github.com/OpenNWP/GAME for more information.\n");
    printf("%s", stars);
	printf("What you want to do:\n");
	printf("Run_id:\t\t\t\t%s\n", RUN_ID);
	printf("Run time span:\t\t\t%d s\n", TOTAL_RUN_SPAN);
	printf("Geo properties file:\t\t%s\n", GEO_PROP_FILE);
	printf("Initialization state file:\t%s\n", INIT_STATE_FILE);
	printf("Start year:\t\t\t%d\n", year);
	printf("Start month:\t\t\t%d\n", month);
	printf("Start day:\t\t\t%d\n", day);
	printf("Start hour:\t\t\t%d\n", hour);
	printf("%s", stars);
	printf("Dynamics configuration:\n");
	printf("Number of layers: %d\n", NO_OF_LAYERS);
	printf("Number of scalar data points per layer: %d\n", NO_OF_SCALARS_H);
	printf("Number of horizontal vectors per layer: %d\n", NO_OF_VECTORS_H);
	printf("Number of scalar data points: %d\n", NO_OF_SCALARS);
	printf("Number of vectors: %d\n", NO_OF_VECTORS);
	printf("Number of data points: %d\n", NO_OF_SCALARS + NO_OF_VECTORS);
	printf("Ratio of the slow to the fast time step: %d\n", config -> slow_fast_ratio);
	if (VERT_GRID_TYPE == 0)
	{
		printf("Terrain handling: terrain following coordinates\n");
		printf("Number of layers following orography: %d\n", grid -> no_of_oro_layers);
	}
	if (VERT_GRID_TYPE == 1)
	{
		printf("Terrain handling: block structure\n");
	}
	if (config -> momentum_diff_h == 0)
	{
		printf("Horizontal momentum diffusion is turned off.\n");
	}
	if (config -> momentum_diff_h == 1)
	{
		printf("Horizontal momentum diffusion is turned on.\n");
	}
	if (config -> momentum_diff_v == 0)
	{
		printf("Vertical momentum diffusion is turned off.\n");
	}
	if (config -> momentum_diff_v == 1)
	{
		printf("Vertical momentum diffusion is turned on.\n");
	}
	if (config -> explicit_boundary_layer == 0)
	{
		printf("Explicit boundary layer friction is turned off.\n");
	}
	if (config -> explicit_boundary_layer == 1)
	{
		printf("Explicit boundary layer friction is turned on.\n");
	}
	if (config -> temperature_diff_h == 0)
	{
		printf("Horizontal temperature diffusion is turned off.\n");
	}
	else
	{
		printf("Horizontal temperature diffusion is turned on.\n");
	}
	if (config -> temperature_diff_v == 0)
	{
		printf("Vertical temperature diffusion is turned off.\n");
	}
	else
	{
		printf("Vertical temperature diffusion is turned on.\n");
	}
	if (config -> tracer_diff_h == 0)
	{
		printf("Horizontal tracer diffusion is turned off.\n");
	}
	else
	{
		printf("Horizontal tracer diffusion is turned on.\n");
	}
	if (config -> tracer_diff_v == 0)
	{
		printf("Vertical tracer diffusion is turned off.\n");
	}
	else
	{
		printf("Vertical tracer diffusion is turned on.\n");
	}
	printf("Horizontal diffusion Smagorinsky factor: %lf.\n", config -> diff_h_smag_fac);
	printf("Background shear: %lf 1/s.\n", config -> shear_bg);
	printf("Swamp layer starts at %lf m.\n", config -> damping_start_height_over_toa*toa);
	printf("Maximum swamp layer damping coefficient: %lf 1/s.\n", config -> damping_coeff_max);
	printf("%s", stars);
	
	printf("Physics configuration:\n");
	printf("Number of constituents: %d\n", NO_OF_CONSTITUENTS);
	printf("Number of condensed constituents: %d\n", NO_OF_CONDENSED_CONSTITUENTS);
	printf("Number of gaseous constituents: %d\n", NO_OF_GASEOUS_CONSTITUENTS);
	if (NO_OF_CONSTITUENTS != 1 && NO_OF_CONSTITUENTS != 6)
	{
		printf("Error: NO_OF_CONSTITUENTS must be either 1 or 6.\n");
		printf("Aborting.\n");
		exit(1);
	}
	if (config -> assume_lte == 0)
	{
		printf("Not Assuming local thermodynamic equilibrium.\n");
	}
	if (config -> assume_lte == 1)
	{
		printf("Assuming local thermodynamic equilibrium.\n");
	}
	if (config -> rad_on == 0)
	{
		printf("Radiation is turned off.\n");
	}
	if (config -> rad_on == 1)
	{
		printf("Radiation is turned on.\n");
	}
	if (config -> rad_on == 2)
	{
		printf("Held-Suarez-forcing is turned on.\n");
	}
	
	printf("%s", stars);
	printf("I/O configuration:\n");
	printf("Output written in intervals of %d s\n", WRITE_OUT_INTERVAL);
	if (config_io -> grib_output_switch == 0)
	{
		printf("Grib output is turned off.\n");
	}
	else
	{
		printf("Grib output is turned on.\n");
	}
	if (config_io -> netcdf_output_switch == 0)
	{
		printf("Netcdf output is turned off.\n");
	}
	else
	{
		printf("Netcdf output is turned on.\n");
	}
	if (config_io -> model_level_output_switch == 0)
	{
		printf("Model level output is turned off.\n");
	}
	else
	{
		printf("Model level output is turned on.\n");
	}
	if (config_io -> surface_output_switch == 0)
	{
		printf("Surface output is turned off.\n");
	}
	else
	{
		printf("Surface output is turned on.\n");
	}
	if (config_io -> pressure_level_output_switch == 0)
	{
		printf("Pressure level output is turned off.\n");
	}
	else
	{
		printf("Pressure level output is turned on.\n");
	}
	printf("%s", stars);
	printf("Model is fully configured now. Starting to read external data.\n");
	printf("%s", stars);
    // Reading and processing user input finished.
    
    // Reading external data.
	printf("Reading grid data ...\n");
    set_grid_properties(grid, dualgrid, GEO_PROP_FILE);
    // If we have radiation turned on, we also need soil.
    if (config -> rad_on == 1)
    {
    	set_sfc_properties(grid, SFC_PROP_FILE);
    }
    printf("Grid loaded successfully.\n");
    printf("%s", stars);
    printf("Reading initial state ...\n");
    State *state_old = calloc(1, sizeof(State));
    set_init_data(INIT_STATE_FILE, state_old, grid);
	printf("Initial state loaded successfully.\n");
	printf("%s", stars);
	
	printf("Calculating time step ...\n");
    // calculating the average horizontal resolution
	double eff_hor_res = 0;
	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		eff_hor_res += grid -> normal_distance[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + i];
	}
	eff_hor_res = eff_hor_res/NO_OF_VECTORS_H;
    // delta_t is the time step
    double delta_t = dt_parameter*eff_hor_res*1e-3;
    
    // calculating the mean area of the cells
	int layer_index, h_index;
	double cell_area_sum = 0;
	for (int i = 0; i < NO_OF_LEVELS*NO_OF_SCALARS_H; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		cell_area_sum += grid -> area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
	}
	grid -> mean_area_cell = cell_area_sum/(NO_OF_LEVELS*NO_OF_SCALARS_H);
    
    // some more checks and info
    if (radiation_delta_t < delta_t)
    {
    	printf("It is radiation_delta_t < delta_t.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
	printf("Time step set. Information on CFL-related quantities:\n");
    printf("Fast modes time step: %lf s\n", delta_t);
    printf("Slow modes time step: %lf s\n", config -> slow_fast_ratio*delta_t);
	
	// finding the minimum horizontal grid distance
	double normal_dist_min_hor = eff_hor_res;
	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		if(grid -> normal_distance[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + i] < normal_dist_min_hor)
		{
			normal_dist_min_hor = grid -> normal_distance[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + i];
		}
	}
	// finding the minimum vertical grid distance
	double normal_dist_min_vert = grid -> z_vector[0]/NO_OF_LAYERS;
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		if(grid -> normal_distance[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER - NO_OF_SCALARS_H + i] < normal_dist_min_vert)
		{
			normal_dist_min_vert = grid -> normal_distance[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER - NO_OF_SCALARS_H + i];
		}
	}
	printf("Effective horizontal resolution: %lf km\n", 1e-3*eff_hor_res);
	printf("Minimum horizontal normal distance: %lf km\n", 1e-3*normal_dist_min_hor);
    double max_speed_hor = 100;
	printf("Horizontal advective Courant numer: %lf\n", delta_t/normal_dist_min_hor*max_speed_hor);
    double max_speed_vert = 0.1;
	printf("Vertical advective Courant numer: %lf\n", delta_t/normal_dist_min_vert*max_speed_vert);
    printf("%s", stars);
    printf("It begins.\n");
    printf("%s", stars);
    
    int min_no_of_10m_wind_avg_steps = 600/delta_t;
    double *wind_h_10m_array = calloc(1, min_no_of_10m_wind_avg_steps*NO_OF_VECTORS_H*sizeof(double));
    double t_write = t_init;
    double delta_z, wind_closest, wind_second_closest, wind_gradient;
    double vector_to_minimize[NO_OF_LAYERS];
    int closest_index, second_closest_index;
	for (int h_index = 0; h_index < NO_OF_VECTORS_H; ++h_index)
	{
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			vector_to_minimize[j] = fabs(grid -> z_vector[NO_OF_SCALARS_H + (NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + h_index] + 10
			- grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + h_index]);
		}
		closest_index = find_min_index(vector_to_minimize, NO_OF_LAYERS);
		second_closest_index = closest_index - 1;
		if (grid -> z_vector[NO_OF_SCALARS_H + closest_index*NO_OF_VECTORS_PER_LAYER + h_index]
		> grid -> z_vector[NO_OF_SCALARS_H + (NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + h_index] + 10
		&& closest_index < NO_OF_LAYERS - 1)
		{
			second_closest_index = closest_index + 1;
		}
		delta_z
		= 0.5*(grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]]
		+ grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]]) + 10
		- grid -> z_vector[NO_OF_SCALARS_H + closest_index*NO_OF_VECTORS_PER_LAYER + h_index];		
		wind_closest = state_old -> wind[NO_OF_SCALARS_H + closest_index*NO_OF_VECTORS_PER_LAYER + h_index];
		wind_second_closest = state_old -> wind[NO_OF_SCALARS_H + second_closest_index*NO_OF_VECTORS_PER_LAYER + h_index];
		// calculating the vertical gradient of the horizontal wind
		wind_gradient = (wind_closest - wind_second_closest)
		/(grid -> z_vector[NO_OF_SCALARS_H + closest_index*NO_OF_VECTORS_PER_LAYER + h_index] - grid -> z_vector[NO_OF_SCALARS_H + second_closest_index*NO_OF_VECTORS_PER_LAYER + h_index]);
		// here, for all output time steps, the initial value is used
		for (int time_step_10_m_wind = 0; time_step_10_m_wind < min_no_of_10m_wind_avg_steps; ++time_step_10_m_wind)
		{
			wind_h_10m_array[time_step_10_m_wind*NO_OF_VECTORS_H + h_index]
			= wind_closest + delta_z*wind_gradient;
    	}
    }
    Diagnostics *diagnostics = calloc(1, sizeof(Diagnostics));
	temperature_diagnostics(state_old, grid, diagnostics);
	inner_product(state_old -> wind, state_old -> wind, diagnostics -> v_squared, grid);
    Forcings *forcings = calloc(1, sizeof(Forcings));
	Soil *soil = calloc(1, sizeof(Soil));
	init_soil(soil, diagnostics);
    // writing out the initial state of the model run
    write_out(state_old, wind_h_10m_array, min_no_of_10m_wind_avg_steps, t_init, t_write, diagnostics, forcings, grid, dualgrid, RUN_ID, config_io, config, soil);
    t_write += WRITE_OUT_INTERVAL;
    printf("Run progress: %f h\n", (t_init - t_init)/SECONDS_PER_HOUR);
    double t_0;
    t_0 = t_init;
    int time_step_counter = 0;
    clock_t first_time, second_time;
    first_time = clock();
    if (write_out_mass_integrals == 1)
    {
		write_out_integral(state_old, time_step_counter, grid, dualgrid, diagnostics, 0);
	}
    if (write_out_rhotheta_integral == 1)
    {
		write_out_integral(state_old, time_step_counter, grid, dualgrid, diagnostics, 1);
	}
    if (write_out_energy_integrals == 1)
    {
		write_out_integral(state_old, time_step_counter, grid, dualgrid, diagnostics, 2);
	}
    config -> rad_update = 1;
    if (config -> rad_on == 1)
    {
    	radiation_init();
    }
	
	/*
	Preparation of the actual integration.
    --------------------------------------
    */
    double t_rad_update = t_0;
    int wind_10_m_step_counter = 0;
    State *state_tendency = calloc(1, sizeof(State));
    Irreversible_quantities *irrev = calloc(1, sizeof(Irreversible_quantities));
    State *state_new = calloc(1, sizeof(State));
    linear_combine_two_states(state_old, state_old, state_new, 1, 0, grid);
    State *state_write = calloc(1, sizeof(State));
    
    /*
    This is the loop over the time steps.
    -------------------------------------
    */
    // This is necessary because at the very first step of the model integration, some things are handled differently in the time stepping.
    config -> totally_first_step_bool = 1;
    // this is to store the speed of the model integration
    double speed;
    while (t_0 < t_init + TOTAL_RUN_SPAN + 300)
    {
    	// copying the old state to the new state
    	linear_combine_two_states(state_new, state_old, state_old, 1, 0, grid);
    	
    	/*
    	Checking if the radiative fluxes need to be updated:
    	----------------------------------------------------
    	*/
        if (t_0 <= t_rad_update && t_0 + delta_t >= t_rad_update)
        {
        	config -> rad_update = 1;
        	t_rad_update += radiation_delta_t;
        }
        else
        {
        	config -> rad_update = 0;
    	}
    	
    	// Time step integration.
    	manage_rkhevi(state_old, state_new, soil, grid, dualgrid, state_tendency, diagnostics, forcings, irrev, config, delta_t, t_0, time_step_counter);
    	// This switch can be set to zero now and remains there.
    	config -> totally_first_step_bool = 0;
		time_step_counter += 1;	
		
		/*
		Writing out integrals over the model domain if requested by the user.
		---------------------------------------------------------------------
		*/
		if (write_out_mass_integrals == 1)
        {
			write_out_integral(state_new, time_step_counter, grid, dualgrid, diagnostics, 0);
    	}
		if (write_out_rhotheta_integral == 1)
        {
			write_out_integral(state_new, time_step_counter, grid, dualgrid, diagnostics, 1);
    	}
		if (write_out_energy_integrals == 1)
        {
			write_out_integral(state_new, time_step_counter, grid, dualgrid, diagnostics, 2);
    	}
    	
    	/*
    	Writing the actual output.
    	--------------------------
    	// interpolating to the output time
    	*/
        if(t_0 + delta_t >= t_write && t_0 <= t_write)
        {
            interpolation_t(state_old, state_new, state_write, t_0, t_0 + delta_t, t_write, grid);
    	}
        
        // 5 minutes before the output time, the wind in the lowest layer needs to be collected for 10 m wind diagnostics.
        if (t_0 >= t_write - 300)
        {
        	if (wind_10_m_step_counter < min_no_of_10m_wind_avg_steps)
        	{
		    	for (int h_index = 0; h_index < NO_OF_VECTORS_H; ++h_index)
       			{
   					for (int j = 0; j < NO_OF_LAYERS; ++j)
					{
						vector_to_minimize[j] = fabs(grid -> z_vector[NO_OF_SCALARS_H + (NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + h_index] + 10
						- grid -> z_vector[NO_OF_SCALARS_H + j*NO_OF_VECTORS_PER_LAYER + h_index]);
					}
					closest_index = find_min_index(vector_to_minimize, NO_OF_LAYERS);
					second_closest_index = closest_index - 1;
					if (grid -> z_vector[NO_OF_SCALARS_H + closest_index*NO_OF_VECTORS_PER_LAYER + h_index]
					> grid -> z_vector[NO_OF_SCALARS_H + (NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + h_index] + 10
					&& closest_index < NO_OF_LAYERS - 1)
					{
						second_closest_index = closest_index + 1;
					}
					// the vertical distance between the desired position and the closest layers
					delta_z
					= 0.5*(grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]]
					+ grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]]) + 10
					- grid -> z_vector[NO_OF_SCALARS_H + closest_index*NO_OF_VECTORS_PER_LAYER + h_index];
					// wind in the closest layer
					wind_closest = state_old -> wind[NO_OF_SCALARS_H + closest_index*NO_OF_VECTORS_PER_LAYER + h_index];
					// wind in the other layer
					wind_second_closest = state_old -> wind[NO_OF_SCALARS_H + second_closest_index*NO_OF_VECTORS_PER_LAYER + h_index];
					// calculating the vertical gradient of the horizontal wind
					wind_gradient = (wind_closest - wind_second_closest)
					// the vertical distance between the two layers used for calculating the gradient
					/(grid -> z_vector[NO_OF_SCALARS_H + closest_index*NO_OF_VECTORS_PER_LAYER + h_index]
					- grid -> z_vector[NO_OF_SCALARS_H + second_closest_index*NO_OF_VECTORS_PER_LAYER + h_index]);
					// interpolating or extrapolating the wind to the desired position
		    		wind_h_10m_array[wind_10_m_step_counter*NO_OF_VECTORS_H + h_index]
		    		= wind_closest + delta_z*wind_gradient;
		    	}
		    	wind_10_m_step_counter += 1;
        	}
        }
        // 5 minutes after the output time, the 10 m wind diagnostics can be executed, so output can actually be written
        if(t_0 + delta_t >= t_write + 300 && t_0 <= t_write + 300)
        {
        	// here, output is actually written
            write_out(state_write, wind_h_10m_array, min_no_of_10m_wind_avg_steps, t_init, t_write, diagnostics, forcings, grid, dualgrid, RUN_ID, config_io, config, soil);
            // setting the next output time
            t_write += WRITE_OUT_INTERVAL;
            
            // Calculating the speed of the model.
            second_time = clock();
        	speed = CLOCKS_PER_SEC*WRITE_OUT_INTERVAL/((double) second_time - first_time);
            printf("Current speed: %lf\n", speed);
            first_time = clock();
            printf("Run progress: %f h\n", (t_0 + delta_t - t_init)/SECONDS_PER_HOUR);
            
            // resetting the wind in the lowest layer to zero
            for (int i = 0; i < min_no_of_10m_wind_avg_steps*NO_OF_VECTORS_H; ++i)
            {
            	wind_h_10m_array[i] = 0;
        	}
            wind_10_m_step_counter = 0;
        }
        // giving the user information on the run progress
        printf("Step %d completed.\n", time_step_counter);
        
    	// updating the model time
        t_0 += delta_t;
    }
    
    /*
    Clean-up.
    ---------
    */
    free(month_string);
    free(day_string);
    free(hour_string);
   	free(RUN_ID);
    free(irrev);
    free(config);
    free(config_io);
    free(diagnostics);
    free(forcings);
    free(state_tendency);
    free(grid);
    free(dualgrid);
    free(state_old);
    free(state_new);
    free(state_write);
    printf("%s", stars);
    free(stars);
    clock_t end = clock();
    speed = CLOCKS_PER_SEC*(TOTAL_RUN_SPAN + 300)/((double) end - begin);
    printf("Average speed: %lf\n", speed);
    printf("GAME over.\n");
    return 0;
}
















