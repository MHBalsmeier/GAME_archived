/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "../spatial_operators/spatial_operators.h"
#include "diagnostics.h"
#include <stdio.h>

int temperature_diagnostics(Scalar_field density_entropy, Scalar_field density, Tracer_densities tracer_densities, Scalar_field temperature)
{
    double exner_pressure, pot_temp, condensates_density_sum, density_d_micro_value, density_v_micro_value;
    int layer_index, h_index;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
    	layer_index = i/NUMBER_OF_SCALARS_H;
    	h_index = i - layer_index*NUMBER_OF_SCALARS_H;
    	condensates_density_sum = calc_condensates_density_sum(layer_index, h_index, tracer_densities);
    	pot_temp = pot_temp_diagnostics_single_value(density_entropy[i], density[i], tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i], condensates_density_sum);
    	density_d_micro_value = calc_micro_density(density[i], condensates_density_sum);
    	density_v_micro_value = calc_micro_density(tracer_densities[NUMBER_OF_CONDENSATED_TRACERS*NUMBER_OF_SCALARS + i], condensates_density_sum);
        exner_pressure = exner_pressure_diagnostics_single_value(density_d_micro_value, density_v_micro_value, pot_temp);
        temperature[i] = temperature_diagnostics_single_value(exner_pressure, pot_temp);
    }
    return 0;
}

double temperature_diagnostics_single_value(double exner_pressure_value, double pot_temp_value)
{
	double result = exner_pressure_value*pot_temp_value;
	return result;
}
