/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

/*
In this file, the model run and I/O configurations can be set, which are not accessible via the run script.
*/

// This is for the Klemp (2008) upper boundary damping layer.
int get_damping_layer_properties(double *damping_start_height_over_toa, double *damping_coeff_max)
{
	// This is where the damping starts in relation to the TOA. 0.75, for example, means that the upper 25 % of the atmosphere are affected by Klemp damping.
	*damping_start_height_over_toa = 0.75;
	// The maximum damping coefficient (the damping coefficient increases towards the TOA).
	*damping_coeff_max = 0.2;
	return 0;
}

// Wether or not horizontal wind divergence shall be written out.
int ask_for_divergence_output(int *write_out_divv_h)
{
	*write_out_divv_h = 1;
	return 0;
}

// This function returns the pressure levels for the synoptic output.
int get_synoptic_pressure_levels(double synoptic_pressure_levels[])
{
	synoptic_pressure_levels[0] = 20000;
	synoptic_pressure_levels[1] = 30000;
	synoptic_pressure_levels[2] = 50000;
	synoptic_pressure_levels[3] = 70000;
	synoptic_pressure_levels[4] = 85000;
	synoptic_pressure_levels[5] = 92500;
	return 0;
}

// This function returns the flight levels for the aviation output.
int get_flight_levels(double flight_levels[])
{
	flight_levels[0] = 100;
	flight_levels[1] = 150;
	flight_levels[2] = 200;
	flight_levels[3] = 250;
	flight_levels[4] = 300;
	flight_levels[5] = 350;
	flight_levels[6] = 400;
	return 0;
}






