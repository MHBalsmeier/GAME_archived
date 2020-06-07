
#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "diagnostics.h"
#include <stdio.h>

int exner_pressure_diagnostics(Scalar_field density_entropy, Scalar_field density, Add_comp_densities add_comp_densities , Scalar_field exner_pressure)
{
	double pot_temp, condensates_density_sum, density_d_micro_value, density_v_micro_value;
	int layer_index, h_index;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
    	layer_index = i/NUMBER_OF_SCALARS_H;
    	h_index = i - layer_index*NUMBER_OF_SCALARS_H;
    	condensates_density_sum = calc_condensates_density_sum(layer_index, h_index, add_comp_densities);
    	pot_temp = pot_temp_diagnostics_single_value(density_entropy[i], density[i], add_comp_densities[NUMBER_OF_COND_ADD_COMPS*NUMBER_OF_SCALARS + i], condensates_density_sum);
    	density_d_micro_value = calc_micro_density(density[i], condensates_density_sum);
    	density_v_micro_value = calc_micro_density(add_comp_densities[NUMBER_OF_COND_ADD_COMPS*NUMBER_OF_SCALARS + i], condensates_density_sum);
        exner_pressure[i] = exner_pressure_diagnostics_single_value(density_d_micro_value, density_v_micro_value, pot_temp);
    }
    return 0;
}

double exner_pressure_diagnostics_single_value(double density_d_micro_value, double density_v_micro_value, double pot_temp)
{
	double R_h = gas_constant_diagnostics(density_d_micro_value, density_v_micro_value);
	double result = pow(R_h*(density_d_micro_value + density_v_micro_value)*pot_temp/P_0, R_D/C_D_V);
    return result;
}
