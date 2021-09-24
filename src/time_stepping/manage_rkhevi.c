/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file manages the RKHEVI time stepping.
*/

#include "../enum_and_typedefs.h"
#include "../settings.h"
#include "../spatial_operators/spatial_operators.h"
#include "time_stepping.h"
#include "../radiation/radiation.h"
#include "../thermodynamics/thermodynamics.h"
#include "../io/io.h"
#include "../soil/soil.h"
#include <geos95.h>
#include <stdlib.h>
#include <stdio.h>

int manage_rkhevi(State *state_old, State *state_new, Soil *soil, Grid *grid, Dualgrid *dualgrid, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irrev, Config_info *config_info, double delta_t, double time_coordinate, int total_step_counter)
{
	// slow terms (momentum advection and diffusion) update switch
	int slow_update_bool = 0;
	// delta_t_small is the time step of the divergent modes integration
	double delta_t_small = delta_t;
	// check if slow terms have to be updated
	if (fmod(total_step_counter, config_info -> adv_sound_ratio) == 0)
	{
		// set the respective update switch to one
		slow_update_bool = 1;
		// delta_t is the large time step for the advection integration
		delta_t = config_info -> adv_sound_ratio*delta_t_small;
    }
       
	// diagnosing the temperature
	temperature_diagnostics(state_old, grid, diagnostics);
    
    // interaction with soil (only useful if real radiation is on)
    if (config_info -> rad_on == 1)
    {
    	soil_interaction(soil, diagnostics, forcings, grid, delta_t);
    }
    
	/*
	Loop over the RK substeps.
	*/
	int vector_index;
	for (int i = 0; i < 2; ++i)
	{
		/*
		state_old remains unchanged the whole time.
		At i == 0, it is state_old == state_new.
		*/
		
		// 1.) explicit component of the momentum equation
		// -----------------------------------------------
		// Update of the pressure gradient.
		if (i == 0)
		{
			manage_pressure_gradient(state_new, grid, dualgrid, diagnostics, forcings, irrev, config_info);
		}
		// Only the horizontal momentum is a forward tendency.
		vector_tendencies_expl(state_new, state_tendency, grid, dualgrid, diagnostics, forcings, irrev, config_info, slow_update_bool, i, delta_t);
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
		if (config_info -> rad_on > 0 && config_info -> rad_update == 1 && i == 0)
		{
			call_radiation(state_old, soil, grid, dualgrid, state_tendency, diagnostics, forcings, irrev, config_info, delta_t, time_coordinate);
		}
		scalar_tendencies_expl(state_old, state_new, state_tendency, soil, grid, delta_t, diagnostics, forcings, irrev, config_info, i);

		// 3.) Vertical sound wave solver.
		// -------------------------------
		three_band_solver_ver_waves(state_old, state_new, state_tendency, diagnostics, config_info, delta_t, grid, i);
		
		// 4.) Solving the implicit component of the generalized density equations for tracers.
		// ------------------------------------------------------------------------------------
		if (NO_OF_CONSTITUENTS > 1)
		{
			three_band_solver_gen_densitites(state_old, state_new, state_tendency, diagnostics, config_info, delta_t, grid);
		}
    }
    
    // saturation adjustment, calculation of latent heating rates
    moisturizer(state_new, delta_t, diagnostics, irrev, config_info, grid);

	// in this case, a large time step has been taken, which we modify into a small step here
    if (slow_update_bool == 1 && config_info -> adv_sound_ratio > 1)
    {
    	linear_combine_two_states(state_old, state_new, state_new, 1 - delta_t_small/delta_t, delta_t_small/delta_t, grid);
    }
    
    return 0;
}






