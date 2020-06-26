/*
This source file is part of the Global Geophysical Modeling Frame (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"

double gas_constant_diagnostics(double density_d_micro_value, double density_v_micro_value)
{
	double result = R_D*(1 - density_v_micro_value/(density_d_micro_value + density_v_micro_value) + density_v_micro_value/(density_d_micro_value + density_v_micro_value)*M_D/M_V);
	return result;
}
