/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/
/*
In this source file algebraic conversions and calculations of thermodynamic quantities of a moist atmosphere are collected.
indices as usual:
d:	dry
v:	water vapour
h:	humid
*/

#include "../enum_and_typedefs.h"
#include "../settings.h"
#include "../spatial_operators/spatial_operators.h"
#include "diagnostics.h"
#include <stdio.h>
#include <stdlib.h>

int pot_temp_diagnostics_dry(State *state, Scalar_field pot_temp)
{
	/*
	This is only needed for the output.
	*/
	double condensates_density_sum, density_d_value, r_d, c_d_p;
	#pragma omp parallel for private (condensates_density_sum, density_d_value, r_d, c_d_p)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	condensates_density_sum = calc_condensates_density_sum(i, state -> mass_densities);
    	density_d_value = calc_micro_density(state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i], condensates_density_sum);
    	r_d = specific_gas_constants(0);
    	c_d_p = spec_heat_capacities_p_gas(0);
    	pot_temp[i] = state -> temperature_gas[i]*pow(P_0/(density_d_value*r_d*state -> temperature_gas[i]), r_d/c_d_p);
	}
    return 0;
}

double spec_heat_cap_diagnostics_v(State *state, int grid_point_index, Config_info *config_info)
{
	double rho_g = 0;
	int no_of_relevant_constituents = 0;
	if (config_info -> assume_lte == 0)
	{
		no_of_relevant_constituents = NO_OF_GASEOUS_CONSTITUENTS;
		rho_g = density_gas(state, grid_point_index);
	}
	if (config_info -> assume_lte == 1)
	{
		no_of_relevant_constituents = 1;
		rho_g = state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + grid_point_index];
	}
	double result = 0;
	for (int i = 0; i < no_of_relevant_constituents; ++i)
	{
		result += state -> mass_densities[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index]/rho_g*spec_heat_capacities_v_gas(i);
	}
	return result;
}

double spec_heat_cap_diagnostics_p(State *state, int grid_point_index, Config_info *config_info)
{
	double rho_g = 0;
	int no_of_relevant_constituents = 0;
	if (config_info -> assume_lte == 0)
	{
		no_of_relevant_constituents = NO_OF_GASEOUS_CONSTITUENTS;
		rho_g = density_gas(state, grid_point_index);
	}
	if (config_info -> assume_lte == 1)
	{
		no_of_relevant_constituents = 1;
		rho_g = state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + grid_point_index];
	}
	double result = 0;
	for (int i = 0; i < no_of_relevant_constituents; ++i)
	{
		result += state -> mass_densities[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index]/rho_g*spec_heat_capacities_p_gas(i);
	}
	return result;
}

double gas_constant_diagnostics(State *state, int grid_point_index, Config_info *config_info)
{
	double rho_g = 0;
	int no_of_relevant_constituents = 0;
	if (config_info -> assume_lte == 0)
	{
		no_of_relevant_constituents = NO_OF_GASEOUS_CONSTITUENTS;
		rho_g = density_gas(state, grid_point_index);
	}
	if (config_info -> assume_lte == 1)
	{
		no_of_relevant_constituents = 1;
		rho_g = state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + grid_point_index];
	}
	double result = 0;
	for (int i = 0; i < no_of_relevant_constituents; ++i)
	{
		result += state -> mass_densities[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index]/rho_g*specific_gas_constants(i);
	}
	return result;
}

double density_total(State *state, int grid_point_index)
{
	double result = 0;
	for (int i = 0; i < NO_OF_CONSTITUENTS; ++i)
	{
		result += state -> mass_densities[i*NO_OF_SCALARS + grid_point_index];
	}
	return result;
}

double density_gas(State *state, int grid_point_index)
{
	double result = 0;
	// < 1 is temporary, normally it should be < NO_OF_GASEOUS_CONSTITUENTS
	for (int i = 0; i < 1; ++i)
	{
		result += state -> mass_densities[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index];
	}
	return result;
}

double calc_micro_density(double density_macro, double condensates_density_sum)
{
	/*
	In a moist atmosphere one needs to distinguish between the densities with respect to the whole volume and the densities with respect to exclusively the gas phase.
	*/
	double result = density_macro/(1 - condensates_density_sum/RHO_WATER);
	if (result < 0)
	{
		printf("Error: microscopic density is negative.\n");
		printf("Aborting.\n");
		exit(1);
	}
	if (isnan(result))
	{
		printf("Error: microscopic density is nan.\n");
		printf("Aborting.\n");
		exit(1);
	}
	return result;
}

double calc_condensates_density_sum(int scalar_gridpoint_index, Mass_densities mass_densities)
{
	/*
	This is only needed for calculating the "micro densities".
	*/
	double result = 0;
	for (int i = 0; i < NO_OF_CONDENSED_CONSTITUENTS; ++i)
	{
		result += mass_densities[i*NO_OF_SCALARS + scalar_gridpoint_index];
	}
	if (result < 0)
	{
		printf("Error: condensates_density_sum is negative.\n");
		printf("Aborting.\n");
		exit(1);
	}
	if (result >= RHO_WATER)
	{
		printf("Error: condensates_density_sum >= RHO_WATER.\n");
		printf("Aborting.\n");
		exit(1);
	}
	return result;
}





