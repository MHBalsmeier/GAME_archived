#include <math.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "diagnostics.h"

void pressure_diagnostics(Scalar_field pot_temp, Scalar_field density, Scalar_field pressure)
{
	extern double R_d, p_0,kappa;
	for (int i = 0; i<=NUMBER_OF_SCALARS-1; ++i)
	{
		pressure[i] = pow(density[i]*R_d*pot_temp[i],kappa)*pow(1/p_0, kappa-1);
	}
}
