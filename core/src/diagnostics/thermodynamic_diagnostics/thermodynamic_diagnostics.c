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

#include "../../enum_and_typedefs.h"
#include "../../spatial_operators/spatial_operators.h"
#include "../diagnostics.h"
#include <stdio.h>
#include <stdlib.h>

int pot_temp_diagnostics_dry(State *state, Scalar_field pot_temp)
{
	/*
	This is only needed for the output.
	*/
	double condensates_density_sum, density_d_value;
	#pragma omp parallel for private (condensates_density_sum, density_d_value)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	condensates_density_sum = calc_condensates_density_sum(i, state -> tracer_densities);
    	density_d_value = calc_micro_density(state -> density_dry[i], condensates_density_sum);
    	pot_temp[i] = state -> temperature_gas[i]*pow(P_0/(density_d_value*R_D*state -> temperature_gas[i]), R_D/C_D_P);
	}
    return 0;
}

int temperature_diagnostics(State *state_old, State *state_new)
{
    double nominator, denominator, entropy_density_gas_0, entropy_density_gas_1, density_gas_0, density_gas_1, delta_density_gas, delta_entropy_density, temperature_0, specific_entropy_gas_0, specific_entropy_gas_1, c_h_v, c_h_p, R_h, density_d_value, density_v_value;
	#pragma omp parallel for private (nominator, denominator, entropy_density_gas_0, entropy_density_gas_1, density_gas_0, density_gas_1, delta_density_gas, delta_entropy_density, temperature_0, specific_entropy_gas_0, specific_entropy_gas_1, c_h_v, c_h_p, R_h, density_d_value, density_v_value)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	// Difference of the mass densities of the gas phase.
    	density_gas_0 = state_old -> density_dry[i] + 0*state_old -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i];
    	density_gas_1 = state_new -> density_dry[i] + 0*state_new -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i];
    	delta_density_gas = density_gas_1 - density_gas_0;
    	
    	entropy_density_gas_0 = state_old -> entropy_density_dry[i] + 0*state_old -> tracer_entropy_densities[i];
    	entropy_density_gas_1 = state_new -> entropy_density_dry[i] + 0*state_new -> tracer_entropy_densities[i];
    	delta_entropy_density = entropy_density_gas_1 - entropy_density_gas_0;
    	
    	// Specific entropies of the gas phase of the two time steps.
    	specific_entropy_gas_0 = entropy_density_gas_0/density_gas_0;
    	specific_entropy_gas_1 = entropy_density_gas_1/density_gas_1;
    	
    	// The temperature of the gas phase of the old time step.
    	temperature_0 = state_old -> temperature_gas[i];
    	
		// Determining the thermodynamic properties of the gas phase.
    	density_d_value = state_old -> density_dry[i];
    	density_v_value = state_old -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i];
    	c_h_v = spec_heat_cap_diagnostics_v(density_d_value, density_v_value);
    	c_h_p = spec_heat_cap_diagnostics_p(density_d_value, density_v_value);
    	R_h = gas_constant_diagnostics(density_d_value, density_v_value);
    	
    	nominator = c_h_v*density_gas_0*temperature_0 + (R_h*temperature_0 - R_h/c_h_p*specific_entropy_gas_0*temperature_0)*delta_density_gas + R_h/c_h_p*temperature_0*delta_entropy_density;
    	denominator = c_h_v*density_gas_0 + c_h_v/c_h_p*specific_entropy_gas_1*delta_density_gas - c_h_v/c_h_p*delta_entropy_density;
    	state_new -> temperature_gas[i] = nominator/denominator;
    }
    return 0;
}

int temperature_diagnostics_explicit(State *state_old, State *state_tendency, Diagnostics *diagnostics, double delta_t)
{
    double nominator, denominator, entropy_density_gas_0, entropy_density_gas_1, density_gas_0, density_gas_1, delta_density_gas, delta_entropy_density, temperature_0, specific_entropy_gas_0, specific_entropy_gas_1, c_h_v, c_h_p, R_h, density_d_value, density_v_value;
	#pragma omp parallel for private (nominator, denominator, entropy_density_gas_0, entropy_density_gas_1, density_gas_0, density_gas_1, delta_density_gas, delta_entropy_density, temperature_0, specific_entropy_gas_0, specific_entropy_gas_1, c_h_v, c_h_p, R_h, density_d_value, density_v_value)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	// Difference of the mass densities of the gas phase.
    	density_gas_0 = state_old -> density_dry[i] + 0*state_old -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i];
    	density_gas_1 = density_gas_0 + delta_t*(state_tendency -> density_dry[i] + 0*state_tendency -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i]);
    	delta_density_gas = density_gas_1 - density_gas_0;
    	
    	entropy_density_gas_0 = state_old -> entropy_density_dry[i] + 0*state_old -> tracer_entropy_densities[i];
    	entropy_density_gas_1 = entropy_density_gas_0 + delta_t*(state_tendency -> entropy_density_dry[i] + 0*state_tendency -> tracer_entropy_densities[i]);
    	delta_entropy_density = entropy_density_gas_1 - entropy_density_gas_0;
    	
    	// Specific entropies of the gas phase of the two time steps.
    	specific_entropy_gas_0 = entropy_density_gas_0/density_gas_0;
    	specific_entropy_gas_1 = entropy_density_gas_1/density_gas_1;
    	
    	// The temperature of the gas phase of the old time step.
    	temperature_0 = state_old -> temperature_gas[i];
    	
		// Determining the thermodynamic properties of the gas phase.
    	density_d_value = state_old -> density_dry[i];
    	density_v_value = state_old -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i];
    	c_h_v = spec_heat_cap_diagnostics_v(density_d_value, density_v_value);
    	c_h_p = spec_heat_cap_diagnostics_p(density_d_value, density_v_value);
    	R_h = gas_constant_diagnostics(density_d_value, density_v_value);
    	
    	nominator = c_h_v*density_gas_0*temperature_0 + (R_h*temperature_0 - R_h/c_h_p*specific_entropy_gas_0*temperature_0)*delta_density_gas + R_h/c_h_p*temperature_0*delta_entropy_density;
    	denominator = c_h_v*density_gas_0 + c_h_v/c_h_p*specific_entropy_gas_1*delta_density_gas - c_h_v/c_h_p*delta_entropy_density;
    	diagnostics -> temperature_gas_explicit[i] = nominator/denominator;
    }
    return 0;
}

double spec_heat_cap_diagnostics_p(double density_d_value, double density_v_value)
{
	double result = (density_d_value*C_D_P + density_v_value*C_V_P)/(density_d_value + density_v_value);
	return result;
}

double spec_heat_cap_diagnostics_v(double density_d_value, double density_v_value)
{
	double result = (density_d_value*C_D_V + density_v_value*C_V_V)/(density_d_value + density_v_value);
	return result;
}

double gas_constant_diagnostics(double density_d_value, double density_v_value)
{
	double result = R_D*(1 - density_v_value/(density_d_value + density_v_value) + density_v_value/(density_d_value + density_v_value)*M_D/M_V);
	return result;
}

double entropy_constant_diagnostics(double density_d_value, double density_v_value)
{
	double result = (density_d_value*ENTROPY_CONSTANT_D + density_v_value*ENTROPY_CONSTANT_V)/(density_d_value + density_v_value);
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

double calc_condensates_density_sum(int scalar_gridpoint_index, Tracer_densities tracer_densities)
{
	/*
	This is only needed for calculating the "micro densities".
	*/
	double result = 0;
	for (int i = 0; i < NO_OF_CONDENSED_TRACERS; ++i)
	{
		result += tracer_densities[i*NO_OF_SCALARS + scalar_gridpoint_index];
	}
	if (result < -NO_OF_CONDENSED_TRACERS*EPSILON_SECURITY)
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





