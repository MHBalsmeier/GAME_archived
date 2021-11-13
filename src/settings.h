/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

// The number of pressure levels for the pressure level output.
#define NO_OF_PRESSURE_LEVELS 6

int get_gas_contituents_ids(int);
double impl_thermo_weight();
double cloud_droplets_velocity();
double precipitation_droplets_velocity();
int get_pressure_levels(double []);
