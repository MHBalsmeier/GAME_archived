/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

/*
In this file, the model run and I/O configurations can be set, which are not accessible via the run script.
*/

#include "atmostracers.h"
#include "settings.h"
#include "enum_and_typedefs.h"

int get_gas_contituents_ids(int);

// This is for the Klemp (2008) upper boundary damping layer.
int get_damping_layer_properties(double *damping_start_height_over_toa, double *damping_coeff_max)
{
	// This is where the damping starts in relation to the TOA. 0.75, for example, means that the upper 25 % of the atmosphere are affected by Klemp damping.
	*damping_start_height_over_toa = 0.53;
	// The maximum damping coefficient (the damping coefficient increases towards the TOA).
	*damping_coeff_max = 0.25;
	return 0;
}

int get_gas_contituents_ids(int gas_constituent_id)
{
	// This defines the constituents of the gas phase.
	int gas_constituent_ids_vector[NO_OF_GASEOUS_CONSTITUENTS];
	for (int i = 0; i < NO_OF_GASEOUS_CONSTITUENTS; ++i)
	{
		gas_constituent_ids_vector[i] = i;
	}
	return gas_constituent_ids_vector[gas_constituent_id];
}

double get_impl_thermo_weight()
{
	double impl_thermo_weight;
	impl_thermo_weight = spec_heat_capacities_v_gas(0)/spec_heat_capacities_p_gas(0);
	return impl_thermo_weight;
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

// the user should not change anything below here
// ----------------------------------------------

double entropy_constants_gas(int gas_constituent_id)
{
	return entropy_constants_gas_lookup(get_gas_contituents_ids(gas_constituent_id));
}

double mean_particle_masses_gas(int gas_constituent_id)
{
	return mean_particle_masses_gas_lookup(get_gas_contituents_ids(gas_constituent_id));
}

double spec_heat_capacities_v_gas(int gas_constituent_id)
{
	return spec_heat_capacities_v_gas_lookup(get_gas_contituents_ids(gas_constituent_id));
}

double spec_heat_capacities_p_gas(int gas_constituent_id)
{
	return spec_heat_capacities_p_gas_lookup(get_gas_contituents_ids(gas_constituent_id));
}

double specific_gas_constants(int gas_constituent_id)
{
	return specific_gas_constants_lookup(get_gas_contituents_ids(gas_constituent_id));
}



