#include "../enum_and_typedefs.h"
#include "../r_operators/r_operators.h"
#include "diagnostics.h"
#include <stdio.h>

int temperature_diagnostics(Scalar_field entropy_density, Scalar_field density, Scalar_field temperature)
{
    double exner_pressure;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
    {
        exner_pressure = pow(R_D*density[i]*THETA_0*exp(entropy_density[i]/density[i])/P_0, R_D/C_V);
        temperature[i] = exner_pressure*THETA_0*exp(entropy_density[i]/density[i]);
    }
    return 0;
}
