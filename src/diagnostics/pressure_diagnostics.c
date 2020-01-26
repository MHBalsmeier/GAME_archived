#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "diagnostics.h"

void pressure_diagnostics(Scalar_field pot_temp, Scalar_field density, Scalar_field pressure)
{
	for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
		pressure[i] = pow(density[i]*R_D*pot_temp[i], KAPPA)*pow(1/P_0, KAPPA - 1);
}
