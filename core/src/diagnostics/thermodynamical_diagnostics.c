/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
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
#include "../spatial_operators/spatial_operators.h"
#include "diagnostics.h"
#include <stdio.h>
#include <stdlib.h>

int pot_temp_diagnostics(Scalar_field temperature, Scalar_field density, Tracer_densities tracer_densities, Scalar_field pot_temp)
{
	/*
	This is only needed for the output.
	*/
	double condensates_density_sum, density_d_micro_value, density_h_micro_value, density_v_micro_value, R_h;
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	condensates_density_sum = calc_condensates_density_sum(i, tracer_densities);
    	density_d_micro_value = calc_micro_density(density[i], condensates_density_sum);
    	density_v_micro_value = calc_micro_density(tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i], condensates_density_sum);
		R_h = gas_constant_diagnostics(density_d_micro_value, density_v_micro_value);
    	density_h_micro_value = density_d_micro_value + density_v_micro_value;
    	pot_temp[i] = temperature[i]*pow(P_0/(density_h_micro_value*R_h*temperature[i]), R_D/C_D_P);
	}
    return 0;
}

int temperature_diagnostics(Scalar_field density_entropy, Scalar_field density, Tracer_densities tracer_densities, Scalar_field temperature)
{
	/*
	This is actually the equation of state of the gas phase in a moist atmosphere.
	*/
    double condensates_density_sum, density_d_micro_value, density_v_micro_value, density_h_micro_value, R_h, c_h_v, entropy_constant_h;
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	condensates_density_sum = calc_condensates_density_sum(i, tracer_densities);
    	density_d_micro_value = calc_micro_density(density[i], condensates_density_sum);
    	density_v_micro_value = calc_micro_density(tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i], condensates_density_sum);
		R_h = gas_constant_diagnostics(density_d_micro_value, density_v_micro_value);
    	density_h_micro_value = density_d_micro_value + density_v_micro_value;
    	c_h_v = spec_heat_cap_diagnostics_v(density_d_micro_value, density_v_micro_value);
    	entropy_constant_h = entropy_constant_diagnostics(density_d_micro_value, density_v_micro_value);
    	temperature[i] = pow(R_h/P_0, R_h/c_h_v)*exp(-entropy_constant_h/c_h_v)*pow(density_h_micro_value, R_h/c_h_v)*exp(density_entropy[i]/(density_h_micro_value*c_h_v));
    }
    return 0;
}

double spec_heat_cap_diagnostics_p(double density_d_micro_value, double density_v_micro_value)
{
	double result = (density_d_micro_value*C_D_P + density_v_micro_value*C_V_P)/(density_d_micro_value + density_v_micro_value);
	return result;
}

double spec_heat_cap_diagnostics_v(double density_d_micro_value, double density_v_micro_value)
{
	double result = (density_d_micro_value*C_D_V + density_v_micro_value*C_V_V)/(density_d_micro_value + density_v_micro_value);
	return result;
}

double gas_constant_diagnostics(double density_d_micro_value, double density_v_micro_value)
{
	double result = R_D*(1 - density_v_micro_value/(density_d_micro_value + density_v_micro_value) + density_v_micro_value/(density_d_micro_value + density_v_micro_value)*M_D/M_V);
	return result;
}

double entropy_constant_diagnostics(double density_d_micro_value, double density_v_micro_value)
{
	double result = (density_d_micro_value*entropy_constant_d + density_v_micro_value*entropy_constant_v)/(density_d_micro_value + density_v_micro_value);
	return result;
}

double calc_micro_density(double density_macro, double condensates_density_sum)
{
	/*
	In a moist atmosphere one needs to distinguish between the densities with respect to the whole volume and the densities with respect to exclusively the gas phase.
	*/
	double result = density_macro/(1 - condensates_density_sum/RHO_WATER);
	if (result < -EPSILON_TRACERS/(1 - condensates_density_sum/RHO_WATER))
	{
		printf("Error: microscopic density negative.\n");
		exit(1);
	}
	if (isnan(result))
	{
		printf("Error: microscopic density is nan.\n");
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
	for (int i = 0; i < NO_OF_CONDENSATED_TRACERS; ++i)
		result += tracer_densities[i*NO_OF_SCALARS + scalar_gridpoint_index];
	if (result < -NO_OF_CONDENSATED_TRACERS*EPSILON_TRACERS)
	{
		printf("Error: condensates_density_sum negative.\n");
		exit(1);
	}
	if (result >= RHO_WATER)
	{
		printf("Error: condensates_density_sum >= RHO_WATER.\n");
		exit(1);
	}
	return result;
}





