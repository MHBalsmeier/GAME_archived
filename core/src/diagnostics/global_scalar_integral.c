/*
This source file is part of the Global Geophysical Modeling Frame (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"

int global_scalar_integrator(Scalar_field density_gen, Grid *grid, double *result)
{
    *result = 0;
    for (int i = 0; i < NUMBER_OF_SCALARS; ++i)
        *result += density_gen[i]*grid -> volume[i];
    return 0;
}
