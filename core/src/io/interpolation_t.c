/*
This source file is part of the General Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../enum_and_typedefs.h"
#include "../spatial_operators/spatial_operators.h"

int interpolation_t(State *state_0, State *state_p1, State *state_write, double t_0, double t_p1, double t_write)
{
    double weight_0, weight_p1;
    weight_p1 = (t_write - t_0)/(t_p1 - t_0);
    weight_0 = 1 - weight_p1;
    linear_combine_two_states(state_0, state_p1, state_write, weight_0, weight_p1);
    return 0;
}
