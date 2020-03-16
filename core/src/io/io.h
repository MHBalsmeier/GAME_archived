int set_grid_properties(Grid *, Dualgrid *, char[]);
int calc_delta_t(double *, Grid *);
int interpolation_t(State *, State *, State *, double, double, double);
int set_init_data(char[], State *, double *);
int write_out(State *, double, double, char[]);
int find_string_length_from_int(int, int *);
int find_hour_from_time_coord(double, int *, int *, int *, int *, int *, int *, int *);
