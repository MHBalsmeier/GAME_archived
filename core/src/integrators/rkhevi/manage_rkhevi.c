/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../../enum_and_typedefs.h"
#include "../../spatial_operators/spatial_operators.h"
#include "../integrators.h"
#include "../../diagnostics/diagnostics.h"
#include <geos95.h>
#include <stdlib.h>
#include <stdio.h>

int manage_rkhevi(State *state_old, State *state_new, Extrapolation_info *extrapolation_info, Grid *grid, Dualgrid *dualgrid, Scalar_field radiation_tendency, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irreversible_quantities, Config_info *config_info, double delta_t, double time_coordinate, int total_step_counter)
{
	// slow terms (momentum advection and diffusion) update switch
	int slow_update_bool = 0;
	if (fmod(total_step_counter, config_info -> adv_sound_ratio) == 0)
	{
		slow_update_bool = 1;
    }
    
	/*
	Loop over the RK substeps.
	*/
	int layer_index, h_index;
	for (int i = 0; i < 2; ++i)
	{
		/*
		state_old remains unchanged the whole time.
		At i == 0, it is state_old == state_new.
		*/
		
		// 1.) Explicit component of the momentum equation.
		// ------------------------------------------------
		forward_tendencies(state_new, state_tendency, grid, dualgrid, diagnostics, forcings, extrapolation_info, irreversible_quantities, config_info, i, slow_update_bool, delta_t);
        // time stepping for the horizontal momentum can be directly executed
        for (int j = 0; j < NO_OF_VECTORS; ++j)
        {
        	layer_index = j/NO_OF_VECTORS_PER_LAYER;
        	h_index = j - layer_index*NO_OF_VECTORS_PER_LAYER;
        	if (h_index >= NO_OF_SCALARS_H)
        	{
        		state_new -> velocity_gas[j] = state_old -> velocity_gas[j] + delta_t/(i + 1)*state_tendency -> velocity_gas[j];
        	}
        }
		// Horizontal velocity can be considered to be updated from now on.
		
		// 2.) Explicit component of the generalized density equations.
		// ------------------------------------------------------------
		// state_new contains densities and velcoties from different time levels (velocity already updated, densities not yet)
		backward_tendencies(state_new, state_tendency, grid, dualgrid, delta_t, radiation_tendency, diagnostics, forcings, irreversible_quantities, config_info, i, time_coordinate);
		// determining the explicit component of the new temperature
		
		// 3.) A pre-conditioned new temperature field, only containing explicit entropy and mass density tendencies (including diabatic forcings).
		// ----------------------------------------------------------------------------------------------------------------------------------------
		temperature_diagnostics_explicit(state_old, state_tendency, diagnostics, config_info, delta_t);

		// 4.) Vertical sound wave solver.
		// -------------------------------
		three_band_solver_ver_sound_waves(state_old, state_new, state_tendency, diagnostics, config_info, delta_t, grid, i);
		// Vertical velocity can be seen as updated from now on.
		
		// 5.) Solving the implicit component of the generalized density equaitons.
		// ------------------------------------------------------------------------
		three_band_solver_gen_densitites(state_old, state_new, state_tendency, diagnostics, config_info, delta_t, grid);
    }
    return 0;
}








