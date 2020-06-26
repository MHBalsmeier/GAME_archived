/*
This source file is part of the Global Geophysical Modeling Frame (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "coordinator.h"

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
    char *stars  = malloc(82*sizeof(char));
    double cfl_margin = strtof(argv[6], NULL);
    int dissipation_on;
    dissipation_on = strtod(argv[7], NULL);
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
    int diffusion_on;
    diffusion_on = strtod(argv[14], NULL);
    double radiation_delta_t;
    radiation_delta_t = strtof(argv[15], NULL);
    int tracers_dynamics_delta_t_ratio;
    tracers_dynamics_delta_t_ratio = strtod(argv[16], NULL);
    if (tracers_dynamics_delta_t_ratio < 1)
    {
    	printf("It is tracers_dynamics_delta_t_ratio < 1.\n");
    	exit(1);
    }
    for (int i = 0; i < 82 - 1; ++i)
        stars[i] = '*';
    stars[81] = '\n';
    printf("%s", stars);
    printf("*\t\t\t\t\t\t\t\t\t\t*\n");
    printf("*\t\t\t\tThis is the GAME\t\t\t\t*\n");
    printf("*\t\t\tGlobal Atmospheric Modeling Framework\t\t\t*\n");
    printf("*\t\t\t\t\t\t\t\t\t\t*\n");
    printf("%s", stars);
    printf("Use only legal if authorized by Max H. Balsmeier.\n");
    printf("What you want to do:\n");
    printf("operator:\t\t\t%s\n", OPERATOR);
    free(OPERATOR);
    printf("run time span:\t\t\t%d s\n", TOTAL_RUN_SPAN);
    printf("output written in intervals of\t%d s\n", WRITE_OUT_INTERVAL);
    printf("geo properties file:\t\t%s\n", GEO_PROP_FILE);
    printf("initialization state file:\t%s\n", INIT_STATE_FILE);
    printf("output directory:\t\t%s\n", OUTPUT_FOLDER);
    printf("%s", stars);
    printf("configuration information:\n");
    printf("number of layers: %d\n", NUMBER_OF_LAYERS);
    printf("number of layers following orography: %d\n", NUMBER_OF_ORO_LAYERS);
    printf("number of scalar data points per layer: %d\n", NUMBER_OF_SCALARS_H);
    double surface = 4*M_PI*pow(RADIUS, 2);
    double points_per_axis = pow(NUMBER_OF_SCALARS_H, 0.5);
    int eff_hor_res_km = 1e-3*pow(surface, 0.5)/points_per_axis;
    printf("effective horizontal resolution: %d km\n", eff_hor_res_km);
    printf("number of horizontal vectors per layer: %d\n", NUMBER_OF_VECTORS_H);
    printf("number of scalar data points: %d\n", NUMBER_OF_SCALARS);
    printf("number of vectors: %d\n", NUMBER_OF_VECTORS);
    printf("number of data points: %d\n", NUMBER_OF_SCALARS + NUMBER_OF_VECTORS);
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
    double t_init;
    set_init_data(INIT_STATE_FILE, state_init, &t_init);
    free(INIT_STATE_FILE);
    int min_no_of_output_steps = 600/delta_t;
    double *wind_h_lowest_layer_array = calloc(1, min_no_of_output_steps*NUMBER_OF_VECTORS_H*sizeof(double));
    double t_write = t_init;
    int h_index, time_step_10_m_wind;
	for (int i = 0; i < min_no_of_output_steps*NUMBER_OF_VECTORS_H; ++i)
	{
		time_step_10_m_wind = i/NUMBER_OF_VECTORS_H;
		h_index = i - time_step_10_m_wind*NUMBER_OF_VECTORS_H;
    	wind_h_lowest_layer_array[time_step_10_m_wind*NUMBER_OF_VECTORS_H + h_index] = state_init -> velocity_gas[NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_PER_LAYER + h_index];
    }
    write_out(state_init, wind_h_lowest_layer_array, min_no_of_output_steps, t_init, t_write, OUTPUT_FOLDER, grid);
    t_write += WRITE_OUT_INTERVAL;
    printf("run progress: %f h\n", (t_init - t_init)/SECONDS_PER_HOUR);
    double t_0;
    t_0 = t_init;
    double t_write_integral = t_init;
    State *state_0 = calloc(1, sizeof(State));
    *state_0 = *state_init;
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
	int retval;
	Scalar_field *radiation_tendency = calloc(1, sizeof(Scalar_field));
    if (rad_on == 1)
    {
		retval = calc_rad_heating(*radiation_tendency, NUMBER_OF_SCALARS);
		if (retval != 0)
		{
			printf("Error in calc_rad_heating called from main, position 0.\n");
			exit(1);
		}
    }
    else
    {
    	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    		(*radiation_tendency)[i] = 0;
    }
	t_write_integral += delta_t;
	int tracers_update = 1;
    int counter = 0;
    double *tracer_mass_source_rates = calloc(NUMBER_OF_TRACERS*NUMBER_OF_SCALARS, sizeof(double));
    double *tracer_heat_source_rates = calloc(NUMBER_OF_TRACERS*NUMBER_OF_SCALARS, sizeof(double));
    State *state_tendency = calloc(1, sizeof(State));
    Vector_field *mass_dry_flux_density = calloc(1, sizeof(Vector_field));
    Scalar_field *mass_dry_flux_density_divergence = calloc(1, sizeof(Scalar_field));
    Scalar_field *temperature = calloc(1, sizeof(Scalar_field));
    Vector_field *entropy_gas_flux_density = calloc(1, sizeof(Vector_field));
    Scalar_field *entropy_gas_flux_density_divergence = calloc(1, sizeof(Scalar_field));
    Scalar_field *temp_diffusion_heating = calloc(1, sizeof(Scalar_field));
    Vector_field *temp_gradient = calloc(1, sizeof(Vector_field));
    Vector_field *friction_acc = calloc(1, sizeof(Vector_field));
    Scalar_field *heating_diss = calloc(1, sizeof(Scalar_field));
    Scalar_field *specific_entropy = calloc(1, sizeof(Scalar_field));
    Vector_field *downgradient_macroscopic_energy = calloc(1, sizeof(Vector_field));
    Vector_field *abs_curl_tend = calloc(1, sizeof(Vector_field));
    Curl_field *rel_curl = calloc(1, sizeof(Curl_field));
    Curl_field *abs_curl = calloc(1, sizeof(Curl_field));
    Vector_field *pressure_gradient_acc = calloc(1, sizeof(Vector_field));
    Vector_field *specific_entropy_gradient = calloc(1, sizeof(Vector_field));
    Scalar_field *c_h_p_field = calloc(1, sizeof(Scalar_field));
    Scalar_field *macroscopic_energy = calloc(1, sizeof(Scalar_field));
    Scalar_field *pressure_gradient_decel_factor = calloc(1, sizeof(Scalar_field));
    Scalar_field *pressure_gradient_acc_1 = calloc(1, sizeof(Vector_field));
    Scalar_field *diffusion_coeff_numerical_h = calloc(1, sizeof(Scalar_field));
    Scalar_field *diffusion_coeff_numerical_v = calloc(1, sizeof(Scalar_field));
    Vector_field *mass_dry_diffusion_flux_density = calloc(1, sizeof(Vector_field));
    Scalar_field *mass_dry_diffusion_source_rate = calloc(1, sizeof(Scalar_field));
    Vector_field *temperature_flux_density = calloc(1, sizeof(Vector_field));
    Scalar_field *tracer_density = calloc(1, sizeof(Scalar_field));
    Vector_field *tracer_velocity = calloc(1, sizeof(Vector_field));
    Vector_field *tracer_flux_density = calloc(1, sizeof(Vector_field));
    Scalar_field *tracer_flux_density_divergence = calloc(1, sizeof(Scalar_field));
    Scalar_field *tracer_density_temperature = calloc(1, sizeof(Scalar_field));
    Vector_field *tracer_temperature_flux_density = calloc(1, sizeof(Vector_field));
    Scalar_field *tracer_temperature_flux_density_divergence = calloc(1, sizeof(Scalar_field));
    int rad_update = 1;
    manage_time_stepping(state_0, state_p1, delta_t, grid, dualgrid, dissipation_on, rad_update*rad_on, tracers_on, diffusion_on, *radiation_tendency, tracers_update, tracers_dynamics_delta_t_ratio, tracer_mass_source_rates, tracer_heat_source_rates, state_tendency, *mass_dry_flux_density, *mass_dry_flux_density_divergence, *temperature, *entropy_gas_flux_density, *entropy_gas_flux_density_divergence, *temp_diffusion_heating, *temp_gradient, *friction_acc, *heating_diss, *specific_entropy, *rel_curl, *abs_curl, *downgradient_macroscopic_energy, *pressure_gradient_acc, *abs_curl_tend, *specific_entropy_gradient, *c_h_p_field, *macroscopic_energy, *pressure_gradient_decel_factor, *pressure_gradient_acc_1, *diffusion_coeff_numerical_h, *diffusion_coeff_numerical_v, *mass_dry_diffusion_flux_density, *mass_dry_diffusion_source_rate, *temperature_flux_density, *tracer_density, *tracer_velocity, *tracer_flux_density, *tracer_flux_density_divergence, *tracer_density_temperature, *tracer_temperature_flux_density, *tracer_temperature_flux_density_divergence);
    tracers_update = 0;
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
    while (t_0 + delta_t < t_init + TOTAL_RUN_SPAN + 300)
    {
    	if (counter%tracers_dynamics_delta_t_ratio == 0)
    		tracers_update = 1;
    	else
    		tracers_update = 0;
        t_0 += delta_t;
        *state_0 = *state_p1;
        if (t_0 <= t_rad_update && t_0 + delta_t >= t_rad_update)
        {
        	rad_update = 1;
        	t_rad_update += radiation_delta_t;
        }
        else
        	rad_update = 0;
        manage_time_stepping(state_0, state_p1, delta_t, grid, dualgrid, dissipation_on, rad_update*rad_on, tracers_on, diffusion_on, *radiation_tendency, tracers_update, tracers_dynamics_delta_t_ratio, tracer_mass_source_rates, tracer_heat_source_rates, state_tendency, *mass_dry_flux_density, *mass_dry_flux_density_divergence, *temperature, *entropy_gas_flux_density, *entropy_gas_flux_density_divergence, *temp_diffusion_heating, *temp_gradient, *friction_acc, *heating_diss, *specific_entropy, *rel_curl, *abs_curl, *downgradient_macroscopic_energy, *pressure_gradient_acc, *abs_curl_tend, *specific_entropy_gradient, *c_h_p_field, *macroscopic_energy, *pressure_gradient_decel_factor, *pressure_gradient_acc_1, *diffusion_coeff_numerical_h, *diffusion_coeff_numerical_v, *mass_dry_diffusion_flux_density, *mass_dry_diffusion_source_rate, *temperature_flux_density, *tracer_density, *tracer_velocity, *tracer_flux_density, *tracer_flux_density_divergence, *tracer_density_temperature, *tracer_temperature_flux_density, *tracer_temperature_flux_density_divergence);
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
		    	for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
		    		wind_h_lowest_layer_array[wind_10_m_step_counter*NUMBER_OF_VECTORS_H + i] = state_0 -> velocity_gas[NUMBER_OF_VECTORS - NUMBER_OF_VECTORS_PER_LAYER + i];
		    	wind_10_m_step_counter += 1;
        	}
        }
        if(t_0 + delta_t >= t_write + 300 && t_0 <= t_write + 300)
        {
            write_out(state_write, wind_h_lowest_layer_array, min_no_of_output_steps, t_init, t_write, OUTPUT_FOLDER, grid);
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
            for (int i = 0; i < min_no_of_output_steps*NUMBER_OF_VECTORS_H; ++i)
            	wind_h_lowest_layer_array[i] = 0;
            wind_10_m_step_counter = 0;
        }
    	counter += 1;
    }
    free(wind_h_lowest_layer_array);
    free(tracer_temperature_flux_density_divergence);
    free(tracer_temperature_flux_density);
    free(tracer_density_temperature);
    free(tracer_flux_density_divergence);
    free(tracer_flux_density);
    free(tracer_velocity);
    free(tracer_density);
    free(temperature_flux_density);
    free(mass_dry_diffusion_source_rate);
    free(mass_dry_diffusion_flux_density);
    free(diffusion_coeff_numerical_h);
    free(diffusion_coeff_numerical_v);
    free(pressure_gradient_acc_1);
    free(pressure_gradient_decel_factor);
    free(macroscopic_energy);
    free(c_h_p_field);
    free(specific_entropy_gradient);
    free(abs_curl_tend);
    free(pressure_gradient_acc);
    free(downgradient_macroscopic_energy);
    free(abs_curl);
    free(rel_curl);
	free(specific_entropy);
    free(heating_diss);
    free(friction_acc);
    free(temp_gradient);
    free(temp_diffusion_heating);
    free(entropy_gas_flux_density);
    free(temperature);
    free(mass_dry_flux_density_divergence);
    free(mass_dry_flux_density);
    free(entropy_gas_flux_density_divergence);
    free(state_tendency);
    free(tracer_mass_source_rates);
    free(tracer_heat_source_rates);
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
















