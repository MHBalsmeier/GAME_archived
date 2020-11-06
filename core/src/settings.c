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

// thermodynamic quantities
// ------------------------

double entropy_constants_gas(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 2429487178047751925300627872548148580712448.000000; // ((MEAN_MASS_D*exp(5.0/3))/(3*M_PI*H_BAR*H_BAR))
	}
	if (gas_constituent_id == 1)
	{
		result = 1511084890012154487904341578321985168998400.000000; // ((MEAN_MASS_V*exp(5.0/3))/(3*M_PI*H_BAR*H_BAR))
	}
	return result;
}

double mean_particle_masses_gas(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 0.004810e-23;
	}
	if (gas_constituent_id == 1)
	{
		result = 0.002991e-23;
	}
	return result;
}

double spec_heat_capacities_v_gas(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 717.942189;
	}
	if (gas_constituent_id == 1)
	{
		result = 1396.475121;
	}
	return result;
}

double spec_heat_capacities_p_gas(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 1005.0;
	}
	if (gas_constituent_id == 1)
	{
		result = 1858.0;
	}
	return result;
}

double specific_gas_constants(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 287.057811;
	}
	if (gas_constituent_id == 1)
	{
		result = 461.524879;
	}
	return result;
}

// This function returns the pressure levels for the pressure_leveltic output.
int get_pressure_leveltic_pressure_levels(double pressure_leveltic_pressure_levels[])
{
	pressure_leveltic_pressure_levels[0] = 20000;
	pressure_leveltic_pressure_levels[1] = 30000;
	pressure_leveltic_pressure_levels[2] = 50000;
	pressure_leveltic_pressure_levels[3] = 70000;
	pressure_leveltic_pressure_levels[4] = 85000;
	pressure_leveltic_pressure_levels[5] = 92500;
	return 0;
}

// input and output
// ---------------------

// Wether or not horizontal wind divergence shall be written out.
int ask_for_divergence_output(int *write_out_divv_h)
{
	*write_out_divv_h = 1;
	return 0;
}

// This function returns the flight levels for the flight_level output.
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






