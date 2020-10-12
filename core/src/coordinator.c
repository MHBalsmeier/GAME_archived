/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
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
#include "manage_time_stepping/manage_time_stepping.h"
#include "rte-rrtmgp-c.h"
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
    	exit(0);
    }
    free(WRITE_OUT_INTERVAL_PRE);
    int res_id_input;
    res_id_input = strtod(argv[3], NULL);
    if (res_id_input != RES_ID)
    {
    	printf("You demanded res_id = %d in your input file, but the model has been compiled with RES_ID = %d.\n", res_id_input, RES_ID);
    	printf("Recompile with RES_ID = %d or choose an executable which has been compiled with RES_ID = %d.\n", res_id_input, res_id_input);
    	printf("Aborting.\n");
    	exit(0);
    }
    int no_of_layers_input;
	no_of_layers_input = strtod(argv[4], NULL);
    if (no_of_layers_input != NO_OF_LAYERS)
    {
    	printf("You demanded no_of_layers = %d in your input file, but the model has been compiled with NO_OF_LAYERS = %d.\n", no_of_layers_input, NO_OF_LAYERS);
    	printf("Recompile with NO_OF_LAYERS = %d or choose an executable which has been compiled with NO_OF_LAYERS = %d.\n", no_of_layers_input, no_of_layers_input);
    	printf("Aborting.\n");
    	exit(0);
    }
    double cfl_margin = strtof(argv[5], NULL);
    Config_info *config_info = calloc(1, sizeof(Config_info));
    config_info -> momentum_diff = strtod(argv[6], NULL);
    config_info -> rad_on = strtod(argv[7], NULL);
    config_info -> tracers_on = strtod(argv[8], NULL);
    len = strlen(argv[9]);
    char *OPERATOR = malloc((len + 1)*sizeof(char));
    strcpy(OPERATOR, argv[9]);
    int write_out_dry_mass_integral;
    write_out_dry_mass_integral = strtod(argv[10], NULL);
    int write_out_entropy_integral; 
    write_out_entropy_integral = strtod(argv[11], NULL);
    int write_out_energy_integral;
    write_out_energy_integral = strtod(argv[12], NULL);
    config_info -> temperature_diff_h = strtod(argv[13], NULL);
    double radiation_delta_t;
    radiation_delta_t = strtof(argv[14], NULL);
    int year;
    year = strtod(argv[15], NULL);
    int month;
    month = strtod(argv[16], NULL);
    int day;
    day = strtod(argv[17], NULL);
    int hour;
    hour = strtod(argv[18], NULL);
    config_info -> temperature_diff_v = strtod(argv[19], NULL);
    len = strlen(argv[20]);
    char *RUN_ID = malloc((len + 1)*sizeof(char));
    strcpy(RUN_ID, argv[20]);
    int write_out_linearized_entropy_integral;
    write_out_linearized_entropy_integral = strtod(argv[21], NULL);
    int toa;
	toa = strtod(argv[22], NULL);
    int no_of_oro_layers_input;
	no_of_oro_layers_input = strtod(argv[23], NULL);
	config_info -> mass_dry_diff_h = strtod(argv[24], NULL);
	config_info -> mass_dry_diff_v = strtod(argv[25], NULL);
    Io_config *io_config = calloc(1, sizeof(Io_config));
	io_config -> grib_output_switch = strtod(argv[26], NULL);
	io_config -> netcdf_output_switch = strtod(argv[27], NULL);
	io_config -> synop_output_switch = strtod(argv[28], NULL);
	io_config -> aviation_output_switch = strtod(argv[29], NULL);
	if (io_config -> grib_output_switch == 0 && io_config -> netcdf_output_switch == 0)
	{
		printf("Either grib_output_switch or netcdf_output_switch must be set to 1.\n");
    	printf("Aborting.\n");
		exit(0);
	}
    if (no_of_oro_layers_input != NO_OF_ORO_LAYERS)
    {
    	printf("You demanded no_of_oro_layers = %d in your input file, but the model has been compiled with NO_OF_ORO_LAYERS = %d.\n", no_of_oro_layers_input, NO_OF_ORO_LAYERS);
    	printf("Recompile with NO_OF_ORO_LAYERS = %d or choose an executable which has been compiled with NO_OF_ORO_LAYERS = %d.\n", no_of_oro_layers_input, no_of_oro_layers_input);
    	printf("Aborting.\n");
    	exit(0);
    }
    int TEST_ID;
    TEST_ID = strtod(argv[24], NULL);
    // This sets the ORO_ID (orography ID) as a function of the TEST_ID.
    int ORO_ID;	
	if (TEST_ID == -1)
		ORO_ID = 3;
	if (TEST_ID == 0)
		ORO_ID = 0;
	if (TEST_ID == 1)
		ORO_ID = 1;
	if (TEST_ID == 2 || TEST_ID == 3 || TEST_ID == 4 || TEST_ID == 5)
		ORO_ID = 2;
	if (TEST_ID == 6 || TEST_ID == 7)
		ORO_ID = 3;
	// Determining the name of the grid file from the RES_ID, NO_OF_LAYERS and so on.
    char GEO_PROP_FILE_PRE[200];
	sprintf(GEO_PROP_FILE_PRE, "grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NO_OF_LAYERS, toa, ORO_ID, NO_OF_ORO_LAYERS);
    char GEO_PROP_FILE[strlen(GEO_PROP_FILE_PRE) + 1];
    strcpy(GEO_PROP_FILE, GEO_PROP_FILE_PRE);
	// Determining the name of the init state file from the TEST_ID, RES_ID, NO_OF_LAYERS and so on.
    char INIT_STATE_FILE_PRE[200];
	sprintf(INIT_STATE_FILE_PRE, "input/test_%d_B%dL%dT%d_O%d_OL%d_SCVT.nc", TEST_ID, RES_ID, NO_OF_LAYERS, toa, ORO_ID, NO_OF_ORO_LAYERS);
    char INIT_STATE_FILE[strlen(INIT_STATE_FILE_PRE) + 1];
    strcpy(INIT_STATE_FILE, INIT_STATE_FILE_PRE);
    // determining the time stamp of 2000-01-01
    struct tm offset_t, *p_offset;
    time_t offset_time_pre;
    offset_t.tm_year = 2000 - 1900;
    offset_t.tm_mon = 0;
    offset_t.tm_mday = 1;
    offset_t.tm_hour = 0;
    offset_t.tm_min = 0;
    offset_t.tm_sec = 0;
    offset_t.tm_isdst = -1;
    offset_time_pre = mktime(&offset_t);
    // this is necessary for conversion to UTC
    p_offset = gmtime(&offset_time_pre);
    offset_t = *p_offset;
    offset_time_pre = mktime(&offset_t);
    double offset_time = (double) offset_time_pre;
    // determining the time stamp of the initialization
    struct tm init_t, *p_init;
    time_t init_time_pre;
    init_t.tm_year = year - 1900;
    init_t.tm_mon = month - 1;
    init_t.tm_mday = day;
    init_t.tm_hour = hour;
    init_t.tm_min = 0;
    init_t.tm_sec = 0;
    init_t.tm_isdst = -1;
    init_time_pre = mktime(&init_t);
    // this is necessary for conversion to UTC
    p_init = gmtime(&init_time_pre);
    init_t = *p_init;
    init_time_pre = mktime(&init_t);
    // substracting the C time coordinate of 2000-01-01
    double t_init = (double) init_time_pre - offset_time;
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
    printf("Released under the MIT license, visit https://github.com/MHBalsmeier/game for more information.\n");
	printf("What you want to do:\n");
	printf("operator:\t\t\t%s\n", OPERATOR);
	free(OPERATOR);
	printf("run time span:\t\t\t%d s\n", TOTAL_RUN_SPAN);
	printf("output written in intervals of\t%d s\n", WRITE_OUT_INTERVAL);
	printf("geo properties file:\t\t%s\n", GEO_PROP_FILE);
	printf("initialization state file:\t%s\n", INIT_STATE_FILE);
	printf("Start year:\t\t\t%d\n", year);
	printf("Start month:\t\t\t%d\n", month);
	printf("Start day:\t\t\t%d\n", day);
	printf("Start hour:\t\t\t%d\n", hour);
	printf("%s", stars);
	printf("configuration information:\n");
	printf("number of layers: %d\n", NO_OF_LAYERS);
	printf("number of layers following orography: %d\n", NO_OF_ORO_LAYERS);
	printf("number of scalar data points per layer: %d\n", NO_OF_SCALARS_H);
	double surface = 4*M_PI*pow(RADIUS, 2);
	double points_per_axis = pow(NO_OF_SCALARS_H, 0.5);
	int eff_hor_res_km = 1e-3*pow(surface, 0.5)/points_per_axis;
	printf("effective horizontal resolution: %d km\n", eff_hor_res_km);
	printf("number of horizontal vectors per layer: %d\n", NO_OF_VECTORS_H);
	printf("number of scalar data points: %d\n", NO_OF_SCALARS);
	printf("number of vectors: %d\n", NO_OF_VECTORS);
	printf("number of data points: %d\n", NO_OF_SCALARS + NO_OF_VECTORS);
	if (config_info -> rad_on == 0)
	{
		printf("Radiation is turned off.\n");
	}
	else
	{
		printf("Radiation is turned on.\n");
	}
	if (config_info -> tracers_on == 0)
	{
		printf("Moisture is turned off.\n");
	}
	else
	{
		printf("Moisture is turned on.\n");
	}
	if (config_info -> mass_dry_diff_h == 0)
	{
		printf("Horizontal dry mass diffusion is turned off.\n");
	}
	else
	{
		printf("Horizontal dry mass diffusion is turned on.\n");
	}
	if (config_info -> mass_dry_diff_v == 0)
	{
		printf("Vertical dry mass diffusion is turned off.\n");
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
		printf("Horizontal momentum diffusion is turned on.\n");
	}
	printf("reading grid data and checking ... ");
    set_grid_properties(grid, dualgrid, GEO_PROP_FILE);
    double delta_t;
    calc_delta_t(cfl_margin, &delta_t, grid);
    if (radiation_delta_t < delta_t)
    {
    	printf("It is radiation_delta_t < delta_t.\n");
    	exit(1);
    }
    printf("time step: %lf s\n", delta_t);
    printf("%s", stars);
    printf("It begins.\n");
    printf("%s", stars);
    State *state_old = calloc(1, sizeof(State));
    set_init_data(INIT_STATE_FILE, state_old);
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
    write_out(state_old, wind_h_lowest_layer_array, min_no_of_output_steps, t_init, t_write, diagnostics, forcings, grid, dualgrid, RUN_ID, io_config);
    t_write += WRITE_OUT_INTERVAL;
    printf("run progress: %f h\n", (t_init - t_init)/SECONDS_PER_HOUR);
    double t_0;
    t_0 = t_init;
    double t_write_integral = t_init;
    clock_t first_time, second_time;
    first_time = clock();
    State *state_new = calloc(1, sizeof(State));
    if (write_out_dry_mass_integral == 1)
		write_out_integral(state_old, t_write_integral, RUN_ID, grid, dualgrid, 0);
    if (write_out_entropy_integral == 1)
		write_out_integral(state_old, t_write_integral, RUN_ID, grid, dualgrid, 1);
    if (write_out_energy_integral == 1)
		write_out_integral(state_old, t_write_integral, RUN_ID, grid, dualgrid, 2);
    if (write_out_linearized_entropy_integral == 1)
		write_out_integral(state_old, t_write_integral, RUN_ID, grid, dualgrid, 3);
	Scalar_field *radiation_tendency = calloc(1, sizeof(Scalar_field));
    if (config_info -> rad_on == 1)
    {
		calc_rad_heating(*radiation_tendency, NO_OF_SCALARS);
    }
    else
    {
    	for (int i = 0; i < NO_OF_SCALARS; ++i)
    		(*radiation_tendency)[i] = 0;
    }
	t_write_integral += delta_t;
    int counter = 0;
    State *state_tendency = calloc(1, sizeof(State));
    Interpolate_info *interpolation = calloc(1, sizeof(Interpolate_info));
    Diffusion_info *diffusion = calloc(1, sizeof(Diffusion_info));
    config_info -> rad_update = 1;
    linear_combine_two_states(state_old, state_old, state_new, 1, 0);
    config_info -> totally_first_step_bool = 1;
    manage_time_stepping(state_old, state_new, interpolation, grid, dualgrid, *radiation_tendency, state_tendency, diagnostics, forcings, diffusion, config_info, delta_t);
    counter += 1;
    if (write_out_dry_mass_integral == 1)
		write_out_integral(state_new, t_write_integral, RUN_ID, grid, dualgrid, 0);
    if (write_out_entropy_integral == 1)
		write_out_integral(state_new, t_write_integral, RUN_ID, grid, dualgrid, 1);
    if (write_out_energy_integral == 1)
		write_out_integral(state_new, t_write_integral, RUN_ID, grid, dualgrid, 2);
    if (write_out_linearized_entropy_integral == 1)
		write_out_integral(state_old, t_write_integral, RUN_ID, grid, dualgrid, 3);
	t_write_integral += delta_t;
    State *state_write = calloc(1, sizeof(State));
    double speed;
    config_info -> rad_update = 0;
    double t_rad_update = t_0 + radiation_delta_t;
    int wind_10_m_step_counter = 0;
    int second_write_out_bool = 1;
    MPI_Init(&argc, &argv);
    config_info -> totally_first_step_bool = 0;
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
        	config_info -> rad_update = 0;
            manage_time_stepping(state_old, state_new, interpolation, grid, dualgrid, *radiation_tendency, state_tendency, diagnostics, forcings, diffusion, config_info, delta_t);
		if (write_out_dry_mass_integral == 1)
			write_out_integral(state_new, t_write_integral, RUN_ID, grid, dualgrid, 0);
		if (write_out_entropy_integral == 1)
			write_out_integral(state_new, t_write_integral, RUN_ID, grid, dualgrid, 1);
		if (write_out_energy_integral == 1)
			write_out_integral(state_new, t_write_integral, RUN_ID, grid, dualgrid, 2);
		if (write_out_linearized_entropy_integral == 1)
			write_out_integral(state_old, t_write_integral, RUN_ID, grid, dualgrid, 3);
		t_write_integral += delta_t;
        if(t_0 + delta_t >= t_write && t_0 <= t_write)
            interpolation_t(state_old, state_new, state_write, t_0, t_0 + delta_t, t_write);
        if (t_0 >= t_write - 300)
        {
        	if (wind_10_m_step_counter < min_no_of_output_steps)
        	{
		    	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
		    		wind_h_lowest_layer_array[wind_10_m_step_counter*NO_OF_VECTORS_H + i] = state_old -> velocity_gas[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + i];
		    	wind_10_m_step_counter += 1;
        	}
        }
        if(t_0 + delta_t >= t_write + 300 && t_0 <= t_write + 300)
        {
            write_out(state_write, wind_h_lowest_layer_array, min_no_of_output_steps, t_init, t_write, diagnostics, forcings, grid, dualgrid, RUN_ID, io_config);
            t_write += WRITE_OUT_INTERVAL;
            second_time = clock();
            if (second_write_out_bool == 1)
            {
            	speed = (CLOCKS_PER_SEC*(WRITE_OUT_INTERVAL + 300))/((double) second_time - first_time);
            	second_write_out_bool = 0;
        	}
            else
            	speed = CLOCKS_PER_SEC*WRITE_OUT_INTERVAL/((double) second_time - first_time);
            printf("current speed: %lf\n", speed);
            first_time = clock();
            printf("run progress: %f h\n", (t_0 + delta_t - t_init)/SECONDS_PER_HOUR);
            for (int i = 0; i < min_no_of_output_steps*NO_OF_VECTORS_H; ++i)
            	wind_h_lowest_layer_array[i] = 0;
            wind_10_m_step_counter = 0;
        }
    	counter += 1;
    }
    MPI_Finalize();
   	free(RUN_ID);
    free(diffusion);
    free(config_info);
    free(io_config);
    free(diagnostics);
    free(interpolation);
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
















