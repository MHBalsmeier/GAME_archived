/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, the model run and I/O configurations can be set, which are not accessible via the run script.
*/

#include "atmostracers.h"
#include "settings.h"
#include "enum_and_typedefs.h"

double impl_thermo_weight()
{
	return 0.75;
}

double cloud_droplets_velocity()
{
	return 0.01;
}

double precipitation_droplets_velocity()
{
	return 0.3;
}


// input and output
// ---------------------

// This function returns the pressure levels for the pressure_level output.
int get_pressure_levels(double pressure_levels[])
{
	pressure_levels[0] = 20000;
	pressure_levels[1] = 30000;
	pressure_levels[2] = 50000;
	pressure_levels[3] = 70000;
	pressure_levels[4] = 85000;
	pressure_levels[5] = 92500;
	return 0;
}



