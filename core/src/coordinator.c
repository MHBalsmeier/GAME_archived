#include "coordinator.h"

Grid grid;
Dualgrid dualgrid;

int main(int argc, char *argv[])
{
	char run_cfg_file_pre[strlen(argv[1])];
    for (int i = 0; i < strlen(argv[1]); ++i)
		run_cfg_file_pre[i] = argv[1][i];
	printf("Run configuration file: %s \n", run_cfg_file_pre);
    long  RUN_HRS, WRITE_OUT_INTERVAL_MIN;
    char GEO_PROP_FILE[128];
    char INIT_STATE_FILE[128];
    char OUTPUT_FOLDER[128];
    char line[128];
    FILE *run_cfg_file;
    char dump[28];
    run_cfg_file = fopen(run_cfg_file_pre, "r");
	fgets(line, 128, run_cfg_file);
    sscanf(line,"%s %ld",dump, &RUN_HRS);
    fgets(line, 128, run_cfg_file);
    sscanf(line,"%s %ld",dump, &WRITE_OUT_INTERVAL_MIN);
    fgets(line, 128, run_cfg_file);
    sscanf(line,"%s %s", dump, GEO_PROP_FILE);
    fgets(line, 128, run_cfg_file);
    sscanf(line,"%s %s",dump, INIT_STATE_FILE);
    fgets(line, 128, run_cfg_file);
    sscanf(line,"%s %s",dump, OUTPUT_FOLDER);
    const double WRITE_OUT_INTERVAL = 60*WRITE_OUT_INTERVAL_MIN;
    set_grid_properties(GEO_PROP_FILE);
    long total_run_seconds;
    State state_init;
    double min_dist;
    set_init_data(INIT_STATE_FILE, state_init);
    total_run_seconds = RUN_HRS*SECONDS_PER_HOUR;
    double delta_t = calc_delta_t(RES_ID);
    double write_out_interval = WRITE_OUT_INTERVAL_MIN*60;
    double number_of_steps_pre = total_run_seconds/delta_t;
    long number_of_steps = (long) (ceil(number_of_steps_pre));
    double t_init;
    sscanf(argv[1], "%lf", &t_init);
    double t = t_init;
    write_out(state_init, t_init, t, 0, OUTPUT_FOLDER);
    double t_m1, t_0, t_p1;
	double t_write = t + write_out_interval;
	t_0 = t;
    t_p1 = t + delta_t;
	int write_out_index = 1;
    State state_0 = state_init;
	State tendency_0;
	tendency_0 = tendency(state_0);
    State state_p1 = euler_explicit(state_0, tendency_0, delta_t);
	State state_m1, state_write;
    for (long l = 1; l < number_of_steps; l++)
    {
		printf("Run progress: %f h\n", (t_0 - t)/3600);
		t_m1 = t_0;
		state_m1 = state_0;
		t_0 = t_p1;
		state_0 = state_p1;
		t_p1 += delta_t;
		tendency_0 = tendency(state_0);
		state_p1 = leapfrog(state_m1, state_0, tendency_0, delta_t);
        if(t_p1 >= t_write && t_0 <= t_write)
        {
            state_write = interpolation_t(state_0, state_p1, t_0, t_p1, t_write);
            write_out(state_write, t_init, t_write, write_out_index, OUTPUT_FOLDER);
            t_write += write_out_interval;
			++write_out_index;
        }
    }
    return 0;
}
