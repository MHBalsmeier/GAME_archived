/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"

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
