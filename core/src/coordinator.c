/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
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
#include "time_stepping/time_stepping.h"
#include "rte-rrtmgp-c.h"
#include "geos95.h"
#include <mpi.h>
#include <omp.h>

int main(int argc, char *argv[])
{
    clock_t begin;
    begin = clock();
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
    	exit(1);
    }
    free(WRITE_OUT_INTERVAL_PRE);
    len = strlen(argv[3]);
    char *GEO_PROP_FILE = malloc((len + 1)*sizeof(char));
    strcpy(GEO_PROP_FILE, argv[3]);
    len = strlen(argv[4]);
    char *INIT_STATE_FILE = malloc((len + 1)*sizeof(char));
    strcpy(INIT_STATE_FILE, argv[4]);
    len = strlen(argv[5]);
    char *OUTPUT_FOLDER = malloc((len + 1)*sizeof(char));
    strcpy(OUTPUT_FOLDER, argv[5]);
    double cfl_margin = strtof(argv[6], NULL);
    int momentum_diffusion_on;
    momentum_diffusion_on = strtod(argv[7], NULL);
    int rad_on;
    rad_on = strtod(argv[8], NULL);
    int tracers_on;
    tracers_on = strtod(argv[9], NULL);
    len = strlen(argv[10]);
    char *OPERATOR = malloc((len + 1)*sizeof(char));
    strcpy(OPERATOR, argv[10]);
    int write_out_dry_mass_integral;
    write_out_dry_mass_integral = strtod(argv[11], NULL);
    int write_out_entropy_integral; 
    write_out_entropy_integral = strtod(argv[12], NULL);
    int write_out_energy_integral;
    write_out_energy_integral = strtod(argv[13], NULL);
    int scalar_diffusion_on;
    scalar_diffusion_on = strtod(argv[14], NULL);
    double radiation_delta_t;
    radiation_delta_t = strtof(argv[15], NULL);
    int year;
    year = strtod(argv[16], NULL);
    int month;
    month = strtod(argv[17], NULL);
    int day;
    day = strtod(argv[18], NULL);
    int hour;
    hour = strtod(argv[19], NULL);
    double t_init;
    find_time_coord(year, month, day, hour, 0, 0, 0, &t_init);
    char *stars  = malloc(83*sizeof(char));
    for (int i = 0; i < 81; ++i)
        stars[i] = '*';
    stars[81] = '\n';
    stars[82] = '\0';
    printf("%s", stars);
    printf("*\t\t\t\t\t\t\t\t\t\t*\n");
    printf("*\t\t\t\tThis is the GAME\t\t\t\t*\n");
    printf("*\t\t\tGlobal Atmospheric Modeling Framework\t\t\t*\n");
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
	printf("Start year: %d\n", year);
	printf("Start month: %d\n", month);
	printf("Start day: %d\n", day);
	printf("Start hour: %d\n", hour);
	printf("output directory:\t\t%s\n", OUTPUT_FOLDER);
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
	if (rad_on == 0)
	{
		printf("Radiation is turned off.\n");
	}
	else
	{
		printf("Radiation is turned on.\n");
	}
	if (tracers_on == 0)
	{
		printf("Moisture is turned off.\n");
	}
	else
	{
		printf("Moisture is turned on.\n");
	}
	if (scalar_diffusion_on == 0)
	{
		printf("Scalar diffusion is turned off.\n");
	}
	else
	{
		printf("Scalar diffusion is turned on.\n");
	}
	if (momentum_diffusion_on == 0)
	{
		printf("Momentum diffusion is turned off.\n");
	}
	else
	{
		printf("Momentum diffusion is turned on.\n");
	}
	printf("reading grid data and checking ... ");
    set_grid_properties(grid, dualgrid, GEO_PROP_FILE);
    free(GEO_PROP_FILE);
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
    State *state_init = calloc(1, sizeof(State));
    set_init_data(INIT_STATE_FILE, state_init);
    free(INIT_STATE_FILE);
    int min_no_of_output_steps = 600/delta_t;
    double *wind_h_lowest_layer_array = calloc(1, min_no_of_output_steps*NO_OF_VECTORS_H*sizeof(double));
    double t_write = t_init;
    int h_index, time_step_10_m_wind;
	for (int i = 0; i < min_no_of_output_steps*NO_OF_VECTORS_H; ++i)
	{
		time_step_10_m_wind = i/NO_OF_VECTORS_H;
		h_index = i - time_step_10_m_wind*NO_OF_VECTORS_H;
    	wind_h_lowest_layer_array[time_step_10_m_wind*NO_OF_VECTORS_H + h_index] = state_init -> velocity_gas[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + h_index];
    }
    write_out(state_init, wind_h_lowest_layer_array, min_no_of_output_steps, t_init, t_write, OUTPUT_FOLDER, grid, dualgrid);
    t_write += WRITE_OUT_INTERVAL;
    printf("run progress: %f h\n", (t_init - t_init)/SECONDS_PER_HOUR);
    double t_0;
    t_0 = t_init;
    double t_write_integral = t_init;
    State *state_0 = calloc(1, sizeof(State));
    linear_combine_two_states(state_init, state_init, state_0, 1, 0);
    free(state_init);
    clock_t first_time, second_time;
    first_time = clock();
    State *state_p1 = calloc(1, sizeof(State));
    if (write_out_dry_mass_integral == 1)
		write_out_integral(state_0, t_write_integral, OUTPUT_FOLDER, grid, dualgrid, 0);
    if (write_out_entropy_integral == 1)
		write_out_integral(state_0, t_write_integral, OUTPUT_FOLDER, grid, dualgrid, 1);
    if (write_out_energy_integral == 1)
		write_out_integral(state_0, t_write_integral, OUTPUT_FOLDER, grid, dualgrid, 2);
	Scalar_field *radiation_tendency = calloc(1, sizeof(Scalar_field));
    if (rad_on == 1)
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
    Vector_field *mass_dry_flux_density = calloc(1, sizeof(Vector_field));
    Scalar_field *mass_dry_flux_density_divv = calloc(1, sizeof(Scalar_field));
    Scalar_field *temperature = calloc(1, sizeof(Scalar_field));
    Vector_field *t_tilde_flux_density = calloc(1, sizeof(Vector_field));
    Scalar_field *t_tilde_flux_density_divv = calloc(1, sizeof(Scalar_field));
    Vector_field *temp_gradient = calloc(1, sizeof(Vector_field));
    Scalar_field *specific_entropy = calloc(1, sizeof(Scalar_field));
    Vector_field *pot_vort_tend = calloc(1, sizeof(Vector_field));
    Curl_field *pot_vort = calloc(1, sizeof(Curl_field));
    Vector_field *pressure_gradient_acc = calloc(1, sizeof(Vector_field));
    Vector_field *specific_entropy_gradient = calloc(1, sizeof(Vector_field));
    Scalar_field *c_h_p_field = calloc(1, sizeof(Scalar_field));
    Scalar_field *e_kin_h = calloc(1, sizeof(Scalar_field));
    Vector_field *pressure_gradient_acc_1 = calloc(1, sizeof(Vector_field));
    Vector_field *temperature_flux_density = calloc(1, sizeof(Vector_field));
    Vector_field *temp_gradient_times_c_h_p = calloc(1, sizeof(Vector_field));
    Vector_field *pressure_gradient_acc_old = calloc(1, sizeof(Vector_field));
    Vector_field *e_kin_h_grad = calloc(1, sizeof(Vector_field));
    Scalar_field *wind_field_divv = calloc(1, sizeof(Scalar_field));
    Scalar_field *entropy_gas_flux_density_divv = calloc(1, sizeof(Scalar_field));
    Vector_field *entropy_gas_flux_density = calloc(1, sizeof(Vector_field));
    int rad_update = 1;
    linear_combine_two_states(state_0, state_0, state_p1, 1, 0);
    manage_time_stepping(state_0, state_p1, interpolation, delta_t, grid, dualgrid, momentum_diffusion_on, rad_update*rad_on, tracers_on, scalar_diffusion_on, *radiation_tendency, state_tendency, *mass_dry_flux_density, *mass_dry_flux_density_divv, *temperature, *t_tilde_flux_density, *t_tilde_flux_density_divv, *temp_gradient, *specific_entropy, *pot_vort, *pressure_gradient_acc, *pot_vort_tend, *specific_entropy_gradient, *c_h_p_field, *e_kin_h, *pressure_gradient_acc_1, *temperature_flux_density, *temp_gradient_times_c_h_p, *pressure_gradient_acc_old, *e_kin_h_grad, *wind_field_divv, 1, *entropy_gas_flux_density_divv, *entropy_gas_flux_density, diffusion);
    counter += 1;
    if (write_out_dry_mass_integral == 1)
		write_out_integral(state_p1, t_write_integral, OUTPUT_FOLDER, grid, dualgrid, 0);
    if (write_out_entropy_integral == 1)
		write_out_integral(state_p1, t_write_integral, OUTPUT_FOLDER, grid, dualgrid, 1);
    if (write_out_energy_integral == 1)
		write_out_integral(state_p1, t_write_integral, OUTPUT_FOLDER, grid, dualgrid, 2);
	t_write_integral += delta_t;
    State *state_write = calloc(1, sizeof(State));
    double speed;
    rad_update = 0;
    double t_rad_update = t_0 + radiation_delta_t;
    int wind_10_m_step_counter = 0;
    int second_write_out_bool = 1;
    MPI_Init(&argc, &argv);
    while (t_0 + delta_t < t_init + TOTAL_RUN_SPAN + 300)
    {
        t_0 += delta_t;
    	linear_combine_two_states(state_p1, state_p1, state_0, 1, 0);
        if (t_0 <= t_rad_update && t_0 + delta_t >= t_rad_update)
        {
        	rad_update = 1;
        	t_rad_update += radiation_delta_t;
        }
        else
        	rad_update = 0;
            manage_time_stepping(state_0, state_p1, interpolation, delta_t, grid, dualgrid, momentum_diffusion_on, rad_update*rad_on, tracers_on, scalar_diffusion_on, *radiation_tendency, state_tendency, *mass_dry_flux_density, *mass_dry_flux_density_divv, *temperature, *t_tilde_flux_density, *t_tilde_flux_density_divv, *temp_gradient, *specific_entropy, *pot_vort, *pressure_gradient_acc, *pot_vort_tend, *specific_entropy_gradient, *c_h_p_field, *e_kin_h, *pressure_gradient_acc_1, *temperature_flux_density, *temp_gradient_times_c_h_p, *pressure_gradient_acc_old, *e_kin_h_grad, *wind_field_divv, 0, *entropy_gas_flux_density_divv, *entropy_gas_flux_density, diffusion);
		if (write_out_dry_mass_integral == 1)
			write_out_integral(state_p1, t_write_integral, OUTPUT_FOLDER, grid, dualgrid, 0);
		if (write_out_entropy_integral == 1)
			write_out_integral(state_p1, t_write_integral, OUTPUT_FOLDER, grid, dualgrid, 1);
		if (write_out_energy_integral == 1)
			write_out_integral(state_p1, t_write_integral, OUTPUT_FOLDER, grid, dualgrid, 2);
		t_write_integral += delta_t;
        if(t_0 + delta_t >= t_write && t_0 <= t_write)
            interpolation_t(state_0, state_p1, state_write, t_0, t_0 + delta_t, t_write);
        if (t_0 >= t_write - 300)
        {
        	if (wind_10_m_step_counter < min_no_of_output_steps)
        	{
		    	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
		    		wind_h_lowest_layer_array[wind_10_m_step_counter*NO_OF_VECTORS_H + i] = state_0 -> velocity_gas[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + i];
		    	wind_10_m_step_counter += 1;
        	}
        }
        if(t_0 + delta_t >= t_write + 300 && t_0 <= t_write + 300)
        {
            write_out(state_write, wind_h_lowest_layer_array, min_no_of_output_steps, t_init, t_write, OUTPUT_FOLDER, grid, dualgrid);
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
    free(diffusion);
    free(entropy_gas_flux_density_divv);
    free(entropy_gas_flux_density);
    free(interpolation);
    free(wind_field_divv);
    free(e_kin_h_grad);
    free(pressure_gradient_acc_old);
    free(temp_gradient_times_c_h_p);
    free(wind_h_lowest_layer_array);
    free(temperature_flux_density);
    free(e_kin_h);
    free(c_h_p_field);
    free(specific_entropy_gradient);
    free(pot_vort_tend);
    free(pressure_gradient_acc);
    free(pot_vort);
	free(specific_entropy);
    free(temp_gradient);
    free(t_tilde_flux_density);
    free(temperature);
    free(mass_dry_flux_density_divv);
    free(mass_dry_flux_density);
    free(t_tilde_flux_density_divv);
    free(state_tendency);
    free(radiation_tendency);
    free(grid);
    free(dualgrid);
    free(OUTPUT_FOLDER);
    free(state_0);
    free(state_p1);
    free(state_write);
    printf("%s", stars);
    free(stars);
    clock_t end = clock();
    speed = CLOCKS_PER_SEC*(TOTAL_RUN_SPAN + 300)/((double) end - begin);
    printf("average speed: %lf\n", speed);
    printf("GAME over.\n");
    return 0;
}
















