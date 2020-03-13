int set_grid_properties(Grid *, Dualgrid *, char[]);
int interpolation_t(State *, State *, State *, double, double, double);
int set_init_data(char[], State *);
int write_out(State *, double, double, int, char[]);
double calc_delta_t(int);
int find_string_length_from_int(int, int *);
