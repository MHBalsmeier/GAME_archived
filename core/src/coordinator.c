#include "coordinator.h"

int main(int argc, char *argv[])
{
    clock_t begin;
    begin = clock();
    Grid *grid = malloc(sizeof(Grid));
    Dualgrid *dualgrid = malloc(sizeof(Dualgrid));
    size_t len = strlen(argv[1]);
    char *TOTAL_RUN_SPAN_PRE = malloc((len + 1)*sizeof(char));
    strcpy(TOTAL_RUN_SPAN_PRE, argv[1]);
    long TOTAL_RUN_SPAN = strtol(TOTAL_RUN_SPAN_PRE, NULL, 10);
    free(TOTAL_RUN_SPAN_PRE);
    len = strlen(argv[2]);
    char *WRITE_OUT_INTERVAL_PRE = malloc((len + 1)*sizeof(char));
    strcpy(WRITE_OUT_INTERVAL_PRE, argv[2]);
    long WRITE_OUT_INTERVAL = strtol(WRITE_OUT_INTERVAL_PRE, NULL, 10);
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
    short dissipation_on;
    dissipation_on = strtod(argv[7], NULL);
    for (int i = 0; i < 82 - 1; ++i)
        stars[i] = '*';
    stars[81] = '\n';
    printf("%s", stars);
    printf("*\t\t\t\t\t\t\t\t\t\t*\n");
    printf("*\t\t\t\tThis is the GAME\t\t\t\t*\n");
    printf("*\t\t\tGlobal Geophysical Modelling Frame\t\t\t*\n");
    printf("*\t\t\t\t\t\t\t\t\t\t*\n");
    printf("%s", stars);
    printf("Use only legal if authorized by Max Henrik Balsmeier.\n");
    printf("What you want to do:\n");
    printf("run time span:\t\t\t%ld s\n", TOTAL_RUN_SPAN);
    printf("output written in intervals of\t%ld s\n", WRITE_OUT_INTERVAL);
    printf("geo properties file:\t\t%s\n", GEO_PROP_FILE);
    printf("initialization state file:\t%s\n", INIT_STATE_FILE);
    printf("output directory:\t\t%s\n", OUTPUT_FOLDER);
    printf("%s", stars);
    printf("configuration information:\n");
    printf("number of layers: %d\n", NUMBER_OF_LAYERS);
    printf("number of scalar data points per layer: %d\n", NUMBER_OF_SCALARS_H);
    double surface = 4*M_PI*pow(SEMIMAJOR, 2);
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
    printf("time step: %lf s\n", delta_t);
    printf("%s", stars);
    printf("It begins.\n");
    printf("%s", stars);
    State *state_init = malloc(sizeof(State));
    double t_init;
    set_init_data(INIT_STATE_FILE, state_init, &t_init);
    free(INIT_STATE_FILE);
    write_out(state_init, t_init, 0, OUTPUT_FOLDER, grid);
    printf("run progress: %f h\n", (t_init - t_init)/SECONDS_PER_HOUR);
    double t_0;
    double t_write = t_init + WRITE_OUT_INTERVAL;
    t_0 = t_init;
    int write_out_index = 1;
    State *state_0 = malloc(sizeof(State));
    *state_0 = *state_init;
    free(state_init);
    clock_t first_time, second_time;
    first_time = clock();
    State *state_p1 = malloc(sizeof(State));
    runge_kutta_third_order(state_0, state_p1, delta_t, grid, dualgrid, dissipation_on);
    State *state_write = malloc(sizeof(State));
    double speed;
    if(t_0 + delta_t >= t_write && t_0 <= t_write)
    {
        interpolation_t(state_0, state_p1, state_write, t_0, t_0 + delta_t, t_write);
        write_out(state_write, t_init, t_write, OUTPUT_FOLDER, grid);
        t_write += WRITE_OUT_INTERVAL;
        second_time = clock();
        speed = CLOCKS_PER_SEC*WRITE_OUT_INTERVAL/((double) second_time - first_time);
        printf("current speed: %lf\n", speed);
        first_time = clock();
        printf("run progress: %f h\n", (t_0 + delta_t - t_init)/SECONDS_PER_HOUR);
    }
    while (t_0 + delta_t < t_init + TOTAL_RUN_SPAN)
    {
        t_0 += delta_t;
        *state_0 = *state_p1;
        runge_kutta_third_order(state_0, state_p1, delta_t, grid, dualgrid, dissipation_on);
        if(t_0 + delta_t >= t_write && t_0 <= t_write)
        {
            interpolation_t(state_0, state_p1, state_write, t_0, t_0 + delta_t, t_write);
            write_out(state_write, t_init, t_write, OUTPUT_FOLDER, grid);
            t_write += WRITE_OUT_INTERVAL;
            second_time = clock();
            speed = CLOCKS_PER_SEC*WRITE_OUT_INTERVAL/((double) second_time - first_time);
            printf("current speed: %lf\n", speed);
            first_time = clock();
            printf("run progress: %f h\n", (t_0 + delta_t - t_init)/SECONDS_PER_HOUR);
        }
    }
    free(grid);
    free(dualgrid);
    free(OUTPUT_FOLDER);
    free(state_0);
    free(state_p1);
    free(state_write);
    printf("%s", stars);
    free(stars);
    clock_t end = clock();
    speed = CLOCKS_PER_SEC*TOTAL_RUN_SPAN/((double) end - begin);
    printf("average speed: %lf\n", speed);
    printf("GAME over.\n");
    return 0;
}
