#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "diagnostics.h"
#include <stdio.h>

int exner_pressure_diagnostics(Scalar_field entropy_density, Scalar_field density, Scalar_field exner_pressure)
{
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        exner_pressure[i] = pow(R_D*density[i]*THETA_0*exp(entropy_density[i]/density[i])/P_0, R_D/C_V);
    return 0;
}
