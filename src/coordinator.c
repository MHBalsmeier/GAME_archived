#include "coordinator.h"

double R, R_d, c_v, kappa, min_dist;

int main(int argc, char *argv[])
{
	R = N_A*k_B;
	R_d = R/M_d;
	c_v = c_p - R_d;
	kappa = c_p/c_v;
	min_dist = semimajor;
	set_grid_properties();
	char input_filename[strlen(argv[1])];
	for (int i = 0; i< strlen(argv[1]); ++i)
		input_filename[i] = argv[1][i];
	printf("input file: %s \n",input_filename);
    double total_run_seconds, t, t_init;
    total_run_seconds = RUN_HRS*seconds_per_hour;
    t = find_time_coord(2010, 1, 1, 0, 0, 0, 0);
	t_init = t;
    double delta_t = calc_delta_t(min_dist);
    double write_out_interval = WRITE_OUT_INTERVAL_MINS*60;
    int number_of_steps = total_run_seconds/delta_t + 1;
    if(fmod(total_run_seconds,delta_t) > 0)
	{
		++number_of_steps;
	}
    int check_id = 1;
	State state_init;
    state_init = initializer(input_filename);
    write_out(state_init, t_init, t, 0, input_filename);
    double t_m2, t_m1, t_0;
	double t_write = t + write_out_interval;
    t_m2 = t;
	t_m1 = t;
    t_0 = t + delta_t;
	int write_out_index = 1;
    State state_m1 = state_init;
	State tendency_m1;
	tendency_m1 = tendency(state_m1);
    State state_0 = euler_explicit(state_m1, tendency_m1, delta_t);
	State state_m2, state_write;
    for (long l=2; l<=number_of_steps-1; l=l+1)
    {
		printf("run progress: %f h\n", (t_0-t)/3600);
		t_m2 = t_m1;
		state_m2 = state_m1;
		t_m1 = t_0;
		state_m1 = state_0;
		t_0 = t_0 + delta_t;
		tendency_m1 = tendency(state_m1);
		state_0 = leapfrog(state_m2, state_m1, tendency_m1, delta_t);
        if(t_0 >= t_write && t_m1 < t_write)
        {
            state_write = interpolation_t(state_m1,state_0,t_m1,t_0,t_write);
            write_out(state_write, t_init, t_write, write_out_index, input_filename);
            t_write = t_write + write_out_interval;
			++write_out_index;
        }
    }
    return 0;
}