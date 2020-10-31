/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

// The number of pressure_leveltic pressure levels for the pressure_leveltic output.
#define NO_OF_PRESSURE_LEVELS 6

// The number of flight levels for the flight_level output.
#define NO_OF_FLIGHT_LEVELS 7

int get_damping_layer_properties(double *, double *);
int ask_for_divergence_output(int *);
int get_flight_levels(double []);
int get_pressure_leveltic_pressure_levels(double []);
