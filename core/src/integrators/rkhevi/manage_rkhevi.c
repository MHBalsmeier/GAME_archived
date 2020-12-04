/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "../../spatial_operators/spatial_operators.h"
#include "../integrators.h"
#include "../../diagnostics/diagnostics.h"
#include <geos95.h>
#include <stdlib.h>
#include <stdio.h>

int manage_rkhevi(State *state_old, State *state_new, Interpolation_info *interpolation_info, Grid *grid, Dualgrid *dualgrid, Scalar_field radiation_tendency, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irreversible_quantities, Config_info *config_info, double delta_t, double time_coordinate)
{
    int max_index = find_max_index(state_old -> velocity_gas, NO_OF_VECTORS);
    int min_index = find_min_index(state_old -> velocity_gas, NO_OF_VECTORS);
    double velocity_max = fabs(state_old -> velocity_gas[max_index]);
    if (fabs(state_old -> velocity_gas[min_index]) > velocity_max)
    {
    	velocity_max = fabs(state_old -> velocity_gas[min_index]);
    }
    printf("Maximum velocity: %lf\n", velocity_max);
	/*
	Loop over the RK substeps.
	*/
	int layer_index, h_index;
	double delta_t_rk;
	for (int i = 0; i < config_info -> rk_order; ++i)
	{
		/*
		general remarks:
		-------------------------------------------------------------------------------
		At i == 0, it is state_new == state_old.
		state_old remains unchanged the whole time.
		*/
		
		// 1.) setting the time step of the RK substep
		delta_t_rk = delta_t/(config_info -> rk_order - i);
		
		// 2.) Explicit component of the momentum equation.
		// ----------------------------------------------------------------------------
		forward_tendencies(state_new, state_tendency, grid, dualgrid, diagnostics, forcings, interpolation_info, irreversible_quantities, config_info, i);
        // time stepping for the horizontal momentum can be directly executed
        for (int j = 0; j < NO_OF_VECTORS; ++j)
        {
        	layer_index = j/NO_OF_VECTORS_PER_LAYER;
        	h_index = j - layer_index*NO_OF_VECTORS_PER_LAYER;
        	if (h_index >= NO_OF_SCALARS_H)
        	{
        		state_new -> velocity_gas[j] = state_old -> velocity_gas[j] + delta_t_rk*state_tendency -> velocity_gas[j];
        	}
        }
		// Horizontal velocity can be considered to be updated from now on.
		
		// 3.) Explicit component of the generalized density equations.
		// ----------------------------------------------------------------------------
		backward_tendencies(state_new, interpolation_info, state_tendency, grid, dualgrid, delta_t_rk, radiation_tendency, diagnostics, forcings, irreversible_quantities, config_info, i, time_coordinate);
		// determining the explicit component of the new temperature
		
		// 4.) A pre-conditioned new temperature field, only containing explicit entropy and mass density tendencies (including diabatic forcings).
		// ----------------------------------------------------------------------------
		temperature_diagnostics_explicit(state_old, state_tendency, diagnostics, delta_t_rk);

		// 5.) Vertical sound wave solver.
		// ----------------------------------------------------------------------------
		three_band_solver_ver_sound_waves(state_old, state_tendency, state_new, diagnostics, delta_t_rk, grid);
		// Vertical velocity can be seen as updated from now on.
		
		// 6.) Solving the implicit component of the generalized density equaitons.
		// ----------------------------------------------------------------------------
		three_band_solver_gen_densitites(state_old, state_new, state_tendency, diagnostics, delta_t_rk, grid);
    }
    return 0;
}








