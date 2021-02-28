/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
The main organizes the model, manages the time stepping, calls model output, collects the lowest model layer wind for 10 m wind mean and so on. All the memory needed for integration is allocated and freed here for efficiency (I wonder if this is really relevant).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "enum_and_typedefs.h"
#include "io/io.h"
#include "spatial_operators/spatial_operators.h"
#include "diagnostics/diagnostics.h"
#include "integrators/integrators.h"
#include <mpi.h>

int main(int argc, char *argv[])
{
    clock_t begin = clock();
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
    if (WRITE_OUT_INTERVAL < 900)
    {
    	printf("It is WRITE_OUT_INTERVAL < 900.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
    free(WRITE_OUT_INTERVAL_PRE);
    double cfl_margin = strtof(argv[3], NULL);
    Config_info *config_info = calloc(1, sizeof(Config_info));
    config_info -> momentum_diff = strtod(argv[4], NULL);
    config_info -> rad_on = strtod(argv[5], NULL);
    len = strlen(argv[6]);
    char *OPERATOR = malloc((len + 1)*sizeof(char));
    strcpy(OPERATOR, argv[6]);
    int write_out_dry_mass_integral;
    write_out_dry_mass_integral = strtod(argv[7], NULL);
    int write_out_entropy_integral; 
    write_out_entropy_integral = strtod(argv[8], NULL);
    int write_out_energy_integral;
    write_out_energy_integral = strtod(argv[9], NULL);
    config_info -> temperature_diff_h = strtod(argv[10], NULL);
    double radiation_delta_t;
    radiation_delta_t = strtof(argv[11], NULL);
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
    config_info -> temperature_diff_v = strtod(argv[16], NULL);
    len = strlen(argv[17]);
    char *RUN_ID = malloc((len + 1)*sizeof(char));
    strcpy(RUN_ID, argv[17]);
    int write_out_linearized_entropy_integral;
    write_out_linearized_entropy_integral = strtod(argv[18], NULL);
    int toa;
	toa = strtod(argv[19], NULL);
	int ORO_ID;
	ORO_ID = strtod(argv[20], NULL);
    int IDEAL_INPUT_ID;
    IDEAL_INPUT_ID = strtod(argv[21], NULL);
	config_info -> mass_diff_h = strtod(argv[22], NULL);
	config_info -> mass_diff_v = strtod(argv[23], NULL);
    Io_config *io_config = calloc(1, sizeof(Io_config));
	io_config -> grib_output_switch = strtod(argv[24], NULL);
	io_config -> netcdf_output_switch = strtod(argv[25], NULL);
	io_config -> pressure_level_output_switch = strtod(argv[26], NULL);
	io_config -> flight_level_output_switch = strtod(argv[27], NULL);
	io_config -> model_level_output_switch = strtod(argv[28], NULL);
	io_config -> surface_output_switch = strtod(argv[29], NULL);
	grid -> no_of_oro_layers = strtod(argv[30], NULL);
	int VERT_GRID_TYPE = strtod(argv[31], NULL);
	config_info -> rk_order = strtod(argv[32], NULL);
	config_info -> assume_lte = strtod(argv[33], NULL);
	config_info -> adv_sound_ratio = strtod(argv[34], NULL);
	config_info -> delta_t_between_analyses = strtod(argv[35], NULL);
	if (config_info -> rk_order < 1)
	{
		printf("The Runge-Kutta order must be at least one.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config_info -> assume_lte != 0 && config_info -> assume_lte != 1)
	{
		printf("simplified_moisture_switch must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	// in the case of block-shaped mountains, no lowers follow the orography
	if (VERT_GRID_TYPE == 1)
	{
		grid -> no_of_oro_layers = 0;
	}
	if (io_config -> grib_output_switch == 0 && io_config -> netcdf_output_switch == 0)
	{
		printf("Either grib_output_switch or netcdf_output_switch must be set to 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
    // This sets the ORO_ID (orography ID) as a function of the IDEAL_INPUT_ID.
	if (IDEAL_INPUT_ID == 0 || IDEAL_INPUT_ID == 8 || IDEAL_INPUT_ID == 9)
    {
		ORO_ID = 0;
    }
	if (IDEAL_INPUT_ID == 1)
    {
		ORO_ID = 1;
    }
	if (IDEAL_INPUT_ID == 2 || IDEAL_INPUT_ID == 3 || IDEAL_INPUT_ID == 4 || IDEAL_INPUT_ID == 5)
    {
		ORO_ID = 2;
    }
	if (IDEAL_INPUT_ID == 6 || IDEAL_INPUT_ID == 7 || IDEAL_INPUT_ID == 10 || IDEAL_INPUT_ID == 11 || IDEAL_INPUT_ID == 12)
    {
		ORO_ID = 3;
    }
	// Determining the name of the grid file from the RES_ID, NO_OF_LAYERS and so on.
    char GEO_PROP_FILE_PRE[200];
	sprintf(GEO_PROP_FILE_PRE, "grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NO_OF_LAYERS, toa, ORO_ID, grid -> no_of_oro_layers);
    char GEO_PROP_FILE[strlen(GEO_PROP_FILE_PRE) + 1];
    strcpy(GEO_PROP_FILE, GEO_PROP_FILE_PRE);
	// Determining the name of the init state file from the IDEAL_INPUT_ID, RES_ID, NO_OF_LAYERS and so on.
    char INIT_STATE_FILE_PRE[200];
    // The NWP case.
    if (IDEAL_INPUT_ID == -1)
    {
    	config_info -> nwp_mode = 1;
    	sprintf(INIT_STATE_FILE_PRE, "input/%d%s%s%s_nwp_B%dL%dT%d_O%d_OL%d_SCVT.nc", year, month_string, day_string, hour_string, RES_ID, NO_OF_LAYERS, toa, ORO_ID, grid -> no_of_oro_layers);
    }
    // The idealized input case.
    else
    {
    	config_info -> nwp_mode = 0;
		sprintf(INIT_STATE_FILE_PRE, "input/test_%d_B%dL%dT%d_O%d_OL%d_SCVT.nc", IDEAL_INPUT_ID, RES_ID, NO_OF_LAYERS, toa, ORO_ID, grid -> no_of_oro_layers);
    }
    char INIT_STATE_FILE[strlen(INIT_STATE_FILE_PRE) + 1];
    strcpy(INIT_STATE_FILE, INIT_STATE_FILE_PRE);
    
    /*
    determining the Unix time stamp of the initialization (UTC)
    */
    struct tm init_t;
    init_t.tm_year = year - 1900;
    init_t.tm_mon = month - 1;
    init_t.tm_mday = day;
    init_t.tm_hour = hour;
    init_t.tm_min = 0;
    init_t.tm_sec = 0;
    // turning off DST
    init_t.tm_isdst = 0;
    time_t init_time = mktime(&init_t);
    // converting to double in UTC
    double t_init = (double) init_time + init_t.tm_gmtoff;
    
    // Giving the user some information on the run to about to be executed.
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
    printf("Released under the MIT license, visit https://github.com/AUN4GFD/game for more information.\n");
    printf("%s", stars);
	printf("What you want to do:\n");
	printf("operator:\t\t\t%s\n", OPERATOR);
	printf("run_id:\t\t\t\t%s\n", RUN_ID);
	free(OPERATOR);
	printf("run time span:\t\t\t%d s\n", TOTAL_RUN_SPAN);
	printf("geo properties file:\t\t%s\n", GEO_PROP_FILE);
	printf("initialization state file:\t%s\n", INIT_STATE_FILE);
	printf("Start year:\t\t\t%d\n", year);
	printf("Start month:\t\t\t%d\n", month);
	printf("Start day:\t\t\t%d\n", day);
	printf("Start hour:\t\t\t%d\n", hour);
	printf("%s", stars);
	printf("Dynamics configuration:\n");
	printf("number of layers: %d\n", NO_OF_LAYERS);
	printf("number of scalar data points per layer: %d\n", NO_OF_SCALARS_H);
	printf("number of horizontal vectors per layer: %d\n", NO_OF_VECTORS_H);
	printf("number of scalar data points: %d\n", NO_OF_SCALARS);
	printf("number of vectors: %d\n", NO_OF_VECTORS);
	printf("number of data points: %d\n", NO_OF_SCALARS + NO_OF_VECTORS);
	printf("Runge Kutta order: %d\n", config_info -> rk_order);
	printf("ratio of advective to sound time step: %d\n", config_info -> adv_sound_ratio);
	if (VERT_GRID_TYPE == 0)
	{
		printf("terrain handling: terrain following coordinates\n");
		printf("number of layers following orography: %d\n", grid -> no_of_oro_layers);
	}
	if (VERT_GRID_TYPE == 1)
	{
		printf("terrain handling: block structure\n");
	}
	printf("%s", stars);
	printf("Physics configuration:\n");
	if (config_info -> rad_on == 0)
	{
		printf("Radiation is turned off.\n");
	}
	else
	{
		printf("Radiation is turned on.\n");
	}
	if (config_info -> mass_diff_h == 0)
	{
		printf("Horizontal mass diffusion is turned off.\n");
	}
	else
	{
		printf("Horizontal dry mass diffusion is turned on.\n");
	}
	if (config_info -> mass_diff_v == 0)
	{
		printf("Vertical mass diffusion is turned off.\n");
	}
	else
	{
		printf("Vertical dry mass diffusion is turned on.\n");
	}
	if (config_info -> temperature_diff_h == 0)
	{
		printf("Horizontal temperature diffusion is turned off.\n");
	}
	else
	{
		printf("Horizontal temperature diffusion is turned on.\n");
	}
	if (config_info -> temperature_diff_v == 0)
	{
		printf("Vertical temperature diffusion is turned off.\n");
	}
	else
	{
		printf("Vertical temperature diffusion is turned on.\n");
	}
	if (config_info -> momentum_diff == 0)
	{
		printf("Momentum diffusion is turned off.\n");
	}
	else
	{
		printf("Momentum diffusion is turned on.\n");
	}
	if (config_info -> assume_lte == 0)
	{
		printf("Not Assuming local thermodynamic equilibrium.\n");
	}
	if (config_info -> assume_lte == 1)
	{
		printf("Assuming local thermodynamic equilibrium.\n");
	}
	printf("%s", stars);
	printf("I/O configuration:\n");
	printf("output written in intervals of %d s\n", WRITE_OUT_INTERVAL);
	if (io_config -> grib_output_switch == 0)
	{
		printf("Grib output is turned off.\n");
	}
	else
	{
		printf("Grib output is turned on.\n");
	}
	if (io_config -> netcdf_output_switch == 0)
	{
		printf("Netcdf output is turned off.\n");
	}
	else
	{
		printf("Netcdf output is turned on.\n");
	}
	if (io_config -> model_level_output_switch == 0)
	{
		printf("Model level output is turned off.\n");
	}
	else
	{
		printf("Model level output is turned on.\n");
	}
	if (io_config -> surface_output_switch == 0)
	{
		printf("Surface output is turned off.\n");
	}
	else
	{
		printf("Surface output is turned on.\n");
	}
	if (io_config -> pressure_level_output_switch == 0)
	{
		printf("Pressure level output is turned off.\n");
	}
	else
	{
		printf("Pressure level output is turned on.\n");
	}
	if (io_config -> flight_level_output_switch == 0)
	{
		printf("Flight level output is turned off.\n");
	}
	else
	{
		printf("Flight level output is turned on.\n");
	}
	printf("%s", stars);
	printf("Reading grid data.\n");
    set_grid_properties(grid, dualgrid, GEO_PROP_FILE);
    printf("Grid loaded successfully.\n");
    printf("%s", stars);
    printf("Reading initial state ... ");
    State *state_old = calloc(1, sizeof(State));
    set_init_data(INIT_STATE_FILE, state_old, grid);
    printf("completed.\n");
    // delta_t is the sound time step
    double delta_t;
    calc_delta_t(cfl_margin, &delta_t, grid, state_old, config_info);
    if (radiation_delta_t < delta_t)
    {
    	printf("It is radiation_delta_t < delta_t.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
    // calculating the average horizontal resolution
	double eff_hor_res = 0;
	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		eff_hor_res += grid -> normal_distance[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + i];
	}
	eff_hor_res = eff_hor_res/NO_OF_VECTORS_H;
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
	printf("effective horizontal resolution: %lf km\n", 1e-3*eff_hor_res);
	printf("minimum normal distance: %lf km\n", 1e-3*normal_dist_min_hor);
    printf("sound time step: %lf s\n", delta_t);
    printf("advective time step: %lf s\n", config_info -> adv_sound_ratio*delta_t);
    double max_speed_hor = 100;
	printf("horizontal advective Courant numer: %lf\n", config_info -> adv_sound_ratio*delta_t/normal_dist_min_hor*max_speed_hor);
    double max_speed_vert = 0.1;
	printf("vertical advective Courant numer: %lf\n", config_info -> adv_sound_ratio*delta_t/normal_dist_min_vert*max_speed_vert);
    printf("%s", stars);
    printf("It begins.\n");
    printf("%s", stars);
    // Reading and processing user input finished.
    
    int min_no_of_output_steps = 600/delta_t;
    double *wind_h_lowest_layer_array = calloc(1, min_no_of_output_steps*NO_OF_VECTORS_H*sizeof(double));
    double t_write = t_init;
    int h_index, time_step_10_m_wind;
	for (int i = 0; i < min_no_of_output_steps*NO_OF_VECTORS_H; ++i)
	{
		time_step_10_m_wind = i/NO_OF_VECTORS_H;
		h_index = i - time_step_10_m_wind*NO_OF_VECTORS_H;
    	wind_h_lowest_layer_array[time_step_10_m_wind*NO_OF_VECTORS_H + h_index] = state_old -> velocity_gas[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + h_index];
    }
    Diagnostics *diagnostics = calloc(1, sizeof(Diagnostics));
    Forcings *forcings = calloc(1, sizeof(Forcings));
    write_out(state_old, wind_h_lowest_layer_array, min_no_of_output_steps, t_init, t_write, diagnostics, forcings, grid, dualgrid, RUN_ID, io_config, config_info);
    t_write += WRITE_OUT_INTERVAL;
    printf("run progress: %f h\n", (t_init - t_init)/SECONDS_PER_HOUR);
    double t_0;
    t_0 = t_init;
    int time_step_counter = 0;
    clock_t first_time, second_time;
    first_time = clock();
    if (write_out_dry_mass_integral == 1)
    {
		write_out_integral(state_old, time_step_counter, RUN_ID, grid, dualgrid, diagnostics, 0);
	}
    if (write_out_entropy_integral == 1)
    {
		write_out_integral(state_old, time_step_counter, RUN_ID, grid, dualgrid, diagnostics, 1);
	}
    if (write_out_energy_integral == 1)
    {
		write_out_integral(state_old, time_step_counter, RUN_ID, grid, dualgrid, diagnostics, 2);
	}
    if (write_out_linearized_entropy_integral == 1)
    {
		write_out_integral(state_old, time_step_counter, RUN_ID, grid, dualgrid, diagnostics, 3);
	}
    config_info -> totally_first_step_bool = 1;
    config_info -> rad_update = 1;
	if (config_info -> rad_on == 1)
	{
		radiation_init();
	}
	
	// preparation of the actual integration
    double speed;
    double t_rad_update = t_0;
    int wind_10_m_step_counter = 0;
    int second_write_out_bool = 1;
    MPI_Init(&argc, &argv);
    config_info -> totally_first_step_bool = 0;
    State *state_tendency = calloc(1, sizeof(State));
    Extrapolation_info *extrapolation_info = calloc(1, sizeof(Extrapolation_info));
    Irreversible_quantities *irrev = calloc(1, sizeof(Irreversible_quantities));
    State *state_new = calloc(1, sizeof(State));
    linear_combine_two_states(state_old, state_old, state_new, 1, 0);
	Scalar_field *radiation_tendency = calloc(1, sizeof(Scalar_field));
    State *state_write = calloc(1, sizeof(State));
    
    // this is the loop over the time steps
    while (t_0 + delta_t < t_init + TOTAL_RUN_SPAN + 300)
    {
        t_0 += delta_t;
    	linear_combine_two_states(state_new, state_new, state_old, 1, 0);
        if (t_0 <= t_rad_update && t_0 + delta_t >= t_rad_update)
        {
        	config_info -> rad_update = 1;
        	t_rad_update += radiation_delta_t;
        }
        else
        {
        	config_info -> rad_update = 0;
    	}
    	manage_rkhevi(state_old, state_new, extrapolation_info, grid, dualgrid, *radiation_tendency, state_tendency, diagnostics, forcings, irrev, config_info, delta_t, t_0, time_step_counter);
		time_step_counter += 1;
		if (write_out_dry_mass_integral == 1)
        {
			write_out_integral(state_new, time_step_counter, RUN_ID, grid, dualgrid, diagnostics, 0);
    	}
		if (write_out_entropy_integral == 1)
        {
			write_out_integral(state_new, time_step_counter, RUN_ID, grid, dualgrid, diagnostics, 1);
    	}
		if (write_out_energy_integral == 1)
        {
			write_out_integral(state_new, time_step_counter, RUN_ID, grid, dualgrid, diagnostics, 2);
    	}
		if (write_out_linearized_entropy_integral == 1)
        {
			write_out_integral(state_old, time_step_counter, RUN_ID, grid, dualgrid, diagnostics, 3);
    	}
        if(t_0 + delta_t >= t_write && t_0 <= t_write)
        {
            interpolation_t(state_old, state_new, state_write, t_0, t_0 + delta_t, t_write);
    	}
        if (t_0 >= t_write - 300)
        {
        	if (wind_10_m_step_counter < min_no_of_output_steps)
        	{
		    	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
       			{
		    		wind_h_lowest_layer_array[wind_10_m_step_counter*NO_OF_VECTORS_H + i] = state_old -> velocity_gas[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + i];
		    	}
		    	wind_10_m_step_counter += 1;
        	}
        }
        if(t_0 + delta_t >= t_write + 300 && t_0 <= t_write + 300)
        {
            write_out(state_write, wind_h_lowest_layer_array, min_no_of_output_steps, t_init, t_write, diagnostics, forcings, grid, dualgrid, RUN_ID, io_config, config_info);
            t_write += WRITE_OUT_INTERVAL;
            second_time = clock();
            if (second_write_out_bool == 1)
            {
            	speed = (CLOCKS_PER_SEC*(WRITE_OUT_INTERVAL + 300))/((double) second_time - first_time);
            	second_write_out_bool = 0;
        	}
            else
            {
            	speed = CLOCKS_PER_SEC*WRITE_OUT_INTERVAL/((double) second_time - first_time);
        	}
            printf("current speed: %lf\n", speed);
            first_time = clock();
            printf("run progress: %f h\n", (t_0 + delta_t - t_init)/SECONDS_PER_HOUR);
            for (int i = 0; i < min_no_of_output_steps*NO_OF_VECTORS_H; ++i)
            {
            	wind_h_lowest_layer_array[i] = 0;
        	}
            wind_10_m_step_counter = 0;
        }
        printf("Step %d completed.\n", time_step_counter);
    }
    
    // clean-up
    MPI_Finalize();
    free(month_string);
    free(day_string);
    free(hour_string);
   	free(RUN_ID);
    free(irrev);
    free(config_info);
    free(io_config);
    free(diagnostics);
    free(extrapolation_info);
    free(forcings);
    free(state_tendency);
    free(radiation_tendency);
    free(grid);
    free(dualgrid);
    free(state_old);
    free(state_new);
    free(state_write);
    printf("%s", stars);
    free(stars);
    clock_t end = clock();
    speed = CLOCKS_PER_SEC*(TOTAL_RUN_SPAN + 300)/((double) end - begin);
    printf("average speed: %lf\n", speed);
    printf("GAME over.\n");
    return 0;
}
















