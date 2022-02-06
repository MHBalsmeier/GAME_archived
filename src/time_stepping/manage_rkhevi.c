/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file manages the RKHEVI time stepping.
*/

#include <stdlib.h>
#include <stdio.h>
#include <geos95.h>
#include "../game_types.h"
#include "../spatial_operators/spatial_operators.h"
#include "time_stepping.h"
#include "../radiation/radiation.h"
#include "../thermodynamics/thermodynamics.h"
#include "../io/io.h"

int manage_rkhevi(State *state_old, State *state_new, Grid *grid, Dualgrid *dualgrid, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings,
Irreversible_quantities *irrev, Config *config, double delta_t, double time_coordinate, int total_step_counter)
{
	// slow terms (diffusion) update switch
	int slow_update_bool = 0;
	// check if slow terms have to be updated
	if (fmod(total_step_counter, config -> slow_fast_ratio) == 0)
	{
		// setting the respective update switch to one
		slow_update_bool = 1;
    }
       
	// diagnosing the temperature
	temperature_diagnostics(state_old, grid, diagnostics);
	
	/*
	Loop over the RK substeps.
	*/
	int vector_index;
	for (int rk_step = 0; rk_step < 2; ++rk_step)
	{
		/*
		state_old remains unchanged the whole time.
		At i == 0, it is state_old == state_new.
		*/
		
		// 1.) explicit component of the momentum equation
		// -----------------------------------------------
		// Update of the pressure gradient.
		if (rk_step == 0)
		{
			manage_pressure_gradient(state_new, grid, dualgrid, diagnostics, forcings, irrev, config);
		}
		// Only the horizontal momentum is a forward tendency.
		vector_tendencies_expl(state_new, state_tendency, grid, dualgrid, diagnostics, forcings, irrev, config, slow_update_bool, rk_step, delta_t);
	    // time stepping for the horizontal momentum can be directly executed
	    
	    #pragma omp parallel for private(vector_index)
	    for (int h_index = 0; h_index < NO_OF_VECTORS_H; ++h_index)
	    {
	    	for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
			{
				vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
				state_new -> wind[vector_index] = state_old -> wind[vector_index] + delta_t*state_tendency -> wind[vector_index];
	    	}
	    }
		// Horizontal velocity can be considered to be updated from now on.

		// 2.) Explicit component of the generalized density equations.
		// ------------------------------------------------------------
	    // Radiation is updated here.
		if (config -> rad_on > 0 && config -> rad_update == 1 && rk_step == 0)
		{
			call_radiation(state_old, grid, dualgrid, state_tendency, diagnostics, forcings, irrev, config, delta_t, time_coordinate);
		}
		scalar_tendencies_expl(state_old, state_new, state_tendency, grid, delta_t, diagnostics, forcings, irrev, config, rk_step, slow_update_bool);

		// 3.) Vertical sound wave solver.
		// -------------------------------
		three_band_solver_ver_waves(state_old, state_new, state_tendency, diagnostics, forcings, config, delta_t, grid, rk_step);
		
		// 4.) Solving the implicit component of the generalized density equations for tracers.
		// ------------------------------------------------------------------------------------
		if (NO_OF_CONSTITUENTS > 1)
		{
			three_band_solver_gen_densitites(state_old, state_new, state_tendency, diagnostics, config, delta_t, grid);
		}
		
		// At the second RK step slow terms are never updated.
		slow_update_bool = 0;
    }
    
    // saturation adjustment, calculation of latent heating rates
    moisturizer(state_new, delta_t, diagnostics, irrev, config, grid);
    
    return 0;
}






