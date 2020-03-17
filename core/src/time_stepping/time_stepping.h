int euler_explicit(State *, State *, double, Grid *, Dualgrid *);
int leapfrog(State *, State *, State *, double, Grid *, Dualgrid *);
int adams_bashforth(State *, State *, State *, double, Grid *, Dualgrid *);
int runge_kutta_third_order(State *, State *, double, Grid *, Dualgrid *);
int runge_kutta_fourth_order(State *, State *, double, Grid *, Dualgrid *);
int heun(State *, State *, double, Grid *, Dualgrid *);
