#include "coordinator.h"

int main(int argc, char *argv[])
{
    clock_t begin;
    begin = clock();
    Grid *grid = malloc(sizeof(Grid));
    Dualgrid *dualgrid = malloc(sizeof(Dualgrid));
    char *run_cfg_file_pre = malloc(strlen(argv[1])*sizeof(char));
    for (int i = 0; i < strlen(argv[1]); ++i)
        run_cfg_file_pre[i] = argv[1][i];
    long  RUN_HRS, WRITE_OUT_INTERVAL_MIN;
    char *GEO_PROP_FILE = malloc(128*sizeof(char));
    char *INIT_STATE_FILE = malloc(128*sizeof(char));
    char *OUTPUT_FOLDER = malloc(128*sizeof(char));
    char *line = malloc(128*sizeof(char));
    FILE *run_cfg_file;
    char *dump = malloc(28*sizeof(char));
    run_cfg_file = fopen(run_cfg_file_pre, "r");
    fgets(line, 128, run_cfg_file);
    sscanf(line, "%s %ld",dump, &RUN_HRS);
    fgets(line, 128, run_cfg_file);
    sscanf(line, "%s %ld",dump, &WRITE_OUT_INTERVAL_MIN);
    fgets(line, 128, run_cfg_file);
    sscanf(line, "%s %s", dump, GEO_PROP_FILE);
    fgets(line, 128, run_cfg_file);
    sscanf(line, "%s %s",dump, INIT_STATE_FILE);
    fgets(line, 128, run_cfg_file);
    sscanf(line, "%s %s", dump, OUTPUT_FOLDER);
    free(line);
    double delta_t = calc_delta_t(RES_ID);
    char *stars  = malloc(82*sizeof(char));
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
    printf("run configuration file:\t\t%s\n", run_cfg_file_pre);
    printf("run hours:\t\t\t%ld\n", RUN_HRS);
    printf("output written in intervals of\t%ld min\n", WRITE_OUT_INTERVAL_MIN);
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
    printf("time step: %lf s\n", delta_t);
    printf("%s", stars);
    printf("It begins.\n");
    printf("%s", stars);
    const double WRITE_OUT_INTERVAL = SECONDS_PER_MIN*WRITE_OUT_INTERVAL_MIN;
    set_grid_properties(grid, dualgrid, GEO_PROP_FILE);
    free(GEO_PROP_FILE);
    long total_run_seconds;
    State *state_init = malloc(sizeof(State));
    double min_dist;
    set_init_data(INIT_STATE_FILE, state_init);
    free(INIT_STATE_FILE);
    total_run_seconds = RUN_HRS*SECONDS_PER_HOUR;
    double write_out_interval = WRITE_OUT_INTERVAL_MIN*SECONDS_PER_MIN;
    double t_init;
    sscanf(argv[1], "%lf", &t_init);
    write_out(state_init, t_init, t_init, 0, OUTPUT_FOLDER);
    printf("run progress: %f h\n", (t_init - t_init)/SECONDS_PER_HOUR);
    double t_0, t_p1;
    double t_write = t_init + write_out_interval;
    t_0 = t_init;
    int write_out_index = 1;
    State *state_0 = malloc(sizeof(State));
    *state_0 = *state_init;
    free(state_init);
    State *tendency_0 = malloc(sizeof(State));
    clock_t first_time, second_time;
    first_time = clock();
    tendency(state_0, tendency_0, grid, dualgrid);
    State *state_p1 = malloc(sizeof(State));
    t_p1 = t_init + delta_t/2;
    euler_explicit(state_0, tendency_0, state_p1, delta_t/2);
    State *state_write = malloc(sizeof(State));
    double speed;
    if(t_p1 >= t_write && t_0 <= t_write)
    {
        interpolation_t(state_0, state_p1, state_write, t_0, t_p1, t_write);
        write_out(state_write, t_init, t_write, t_write, OUTPUT_FOLDER);
        t_write += write_out_interval;
        second_time = clock();
        speed = CLOCKS_PER_SEC*WRITE_OUT_INTERVAL/((double) second_time - first_time);
        printf("current speed: %lf\n", speed);
        first_time = clock();
    }
    printf("run progress: %f h\n", (t_p1 - t_init)/SECONDS_PER_HOUR);
    while (t_p1 < t_init + total_run_seconds)
    {
        t_0 = t_p1;
        *state_0 = *state_p1;
        t_p1 += delta_t;
        tendency(state_0, tendency_0, grid, dualgrid);
        euler_explicit(state_0, tendency_0, state_p1, delta_t);
        if(t_p1 >= t_write && t_0 <= t_write)
        {
            interpolation_t(state_0, state_p1, state_write, t_0, t_p1, t_write);
            write_out(state_write, t_init, t_write,t_write, OUTPUT_FOLDER);
            t_write += write_out_interval;
            second_time = clock();
            speed = CLOCKS_PER_SEC*WRITE_OUT_INTERVAL/((double) second_time - first_time);
            printf("current speed: %lf\n", speed);
            first_time = clock();
        }
        printf("run progress: %f h\n", (t_p1 - t_init)/SECONDS_PER_HOUR);
    }
    free(OUTPUT_FOLDER);
    free(run_cfg_file_pre);
    free(state_0);
    free(tendency_0);
    free(state_p1);
    free(state_write);
    free(grid);
    free(dualgrid);
    printf("%s", stars);
    free(stars);
    free(dump);
    clock_t end = clock();
    speed = CLOCKS_PER_SEC*total_run_seconds/((double) end - begin);
    printf("average speed: %lf\n", speed);
    printf("GAME over.\n");
    return 0;
}
