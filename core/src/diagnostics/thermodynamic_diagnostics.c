/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
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

int temperature_diagnostics(State *state_old, State *state_new)
{
    double nominator, denominator, entropy_density_gas_0, entropy_density_gas_1, density_gas_0, density_gas_1, delta_density_gas, delta_entropy_density, temperature_0, specific_entropy_gas_0, specific_entropy_gas_1, c_g_v, c_g_p, R_g;
	#pragma omp parallel for private (nominator, denominator, entropy_density_gas_0, entropy_density_gas_1, density_gas_0, density_gas_1, delta_density_gas, delta_entropy_density, temperature_0, specific_entropy_gas_0, specific_entropy_gas_1, c_g_v, c_g_p, R_g)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	// Difference of the mass densities of the gas phase.
    	for (int j = 0; j < NO_OF_GASEOUS_CONSTITUENTS; ++j)
    	{
			density_gas_0 = state_old -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + j)*NO_OF_SCALARS + i];
			density_gas_1 = state_new -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + j)*NO_OF_SCALARS + i];
		}
    	delta_density_gas = density_gas_1 - density_gas_0;
    	
    	for (int j = 0; j < NO_OF_GASEOUS_CONSTITUENTS; ++j)
    	{
			entropy_density_gas_0 = state_old -> entropy_densities[(NO_OF_CONDENSED_CONSTITUENTS + j)*NO_OF_SCALARS + i];
			entropy_density_gas_1 = state_new -> entropy_densities[(NO_OF_CONDENSED_CONSTITUENTS + j)*NO_OF_SCALARS + i];
		}
    	delta_entropy_density = entropy_density_gas_1 - entropy_density_gas_0;
    	
    	// Specific entropies of the gas phase of the two time steps.
    	specific_entropy_gas_0 = entropy_density_gas_0/density_gas_0;
    	specific_entropy_gas_1 = entropy_density_gas_1/density_gas_1;
    	
    	// The temperature of the gas phase of the old time step.
    	temperature_0 = state_old -> temperature_gas[i];
    	
		// Determining the thermodynamic properties of the gas phase.
    	c_g_v = spec_heat_cap_diagnostics_v(state_old, i);
    	c_g_p = spec_heat_cap_diagnostics_p(state_old, i);
    	R_g = gas_constant_diagnostics(state_old, i);
    	
    	nominator = c_g_v*density_gas_0*temperature_0 + (R_g*temperature_0 - R_g/c_g_p*specific_entropy_gas_0*temperature_0)*delta_density_gas + R_g/c_g_p*temperature_0*delta_entropy_density;
    	denominator = c_g_v*density_gas_0 + c_g_v/c_g_p*specific_entropy_gas_1*delta_density_gas - c_g_v/c_g_p*delta_entropy_density;
    	state_new -> temperature_gas[i] = nominator/denominator;
    }
    return 0;
}

int temperature_diagnostics_explicit(State *state_old, State *state_tendency, Diagnostics *diagnostics, double delta_t)
{
    double nominator, denominator, entropy_density_gas_0, entropy_density_gas_1, density_gas_0, density_gas_1, delta_density_gas, delta_entropy_density, temperature_0, specific_entropy_gas_0, specific_entropy_gas_1, c_g_v, c_g_p, R_g;
	#pragma omp parallel for private (nominator, denominator, entropy_density_gas_0, entropy_density_gas_1, density_gas_0, density_gas_1, delta_density_gas, delta_entropy_density, temperature_0, specific_entropy_gas_0, specific_entropy_gas_1, c_g_v, c_g_p, R_g)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	// Difference of the mass densities of the gas phase.
    	for (int j = 0; j < NO_OF_GASEOUS_CONSTITUENTS; ++j)
    	{
			density_gas_0 = state_old -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + j)*NO_OF_SCALARS + i];
			density_gas_1 = density_gas_0 + delta_t*state_tendency -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + j)*NO_OF_SCALARS + i];
    	}
    	delta_density_gas = density_gas_1 - density_gas_0;
    	
    	for (int j = 0; j < NO_OF_GASEOUS_CONSTITUENTS; ++j)
    	{
			entropy_density_gas_0 = state_old -> entropy_densities[(NO_OF_CONDENSED_CONSTITUENTS + j)*NO_OF_SCALARS + i];
			entropy_density_gas_1 = entropy_density_gas_0 + delta_t*state_tendency -> entropy_densities[(NO_OF_CONDENSED_CONSTITUENTS + j)*NO_OF_SCALARS + i];
    	}
    	delta_entropy_density = entropy_density_gas_1 - entropy_density_gas_0;
    	
    	// Specific entropies of the gas phase of the two time steps.
    	specific_entropy_gas_0 = entropy_density_gas_0/density_gas_0;
    	specific_entropy_gas_1 = entropy_density_gas_1/density_gas_1;
    	
    	// The temperature of the gas phase of the old time step.
    	temperature_0 = state_old -> temperature_gas[i];
    	
		// Determining the thermodynamic properties of the gas phase.
    	c_g_v = spec_heat_cap_diagnostics_v(state_old, i);
    	c_g_p = spec_heat_cap_diagnostics_p(state_old, i);
    	R_g = gas_constant_diagnostics(state_old, i);
    	
    	nominator = c_g_v*density_gas_0*temperature_0 + (R_g*temperature_0 - R_g/c_g_p*specific_entropy_gas_0*temperature_0)*delta_density_gas + R_g/c_g_p*temperature_0*delta_entropy_density;
    	denominator = c_g_v*density_gas_0 + c_g_v/c_g_p*specific_entropy_gas_1*delta_density_gas - c_g_v/c_g_p*delta_entropy_density;
    	diagnostics -> temperature_gas_explicit[i] = nominator/denominator;
    }
    return 0;
}

double spec_heat_cap_diagnostics_v(State *state, int grid_point_index)
{
	double rho_g = density_gas(state, grid_point_index);
	
	double result = 0;
	for (int i = 0; i < NO_OF_GASEOUS_CONSTITUENTS; ++i)
	{
		result += state -> mass_densities[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index]/rho_g*spec_heat_capacities_v_gas(i);
	}
	return result;
}

double spec_heat_cap_diagnostics_p(State *state, int grid_point_index)
{
	double rho_g = density_gas(state, grid_point_index);
	
	double result = 0;
	for (int i = 0; i < NO_OF_GASEOUS_CONSTITUENTS; ++i)
	{
		result += state -> mass_densities[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index]/rho_g*spec_heat_capacities_p_gas(i);
	}
	return result;
}

double gas_constant_diagnostics(State *state, int grid_point_index)
{
	double rho_g = density_gas(state, grid_point_index);
	
	double result = 0;
	for (int i = 0; i < NO_OF_GASEOUS_CONSTITUENTS; ++i)
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
	for (int i = 0; i < NO_OF_GASEOUS_CONSTITUENTS; ++i)
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
	if (result < -EPSILON_SECURITY/(1 - condensates_density_sum/RHO_WATER))
	{
		printf("Error: microscopic density negative.\n");
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
	if (result < -NO_OF_CONDENSED_CONSTITUENTS*EPSILON_SECURITY)
	{
		printf("Error: condensates_density_sum negative.\n");
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





