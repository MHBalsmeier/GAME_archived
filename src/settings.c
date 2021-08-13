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


