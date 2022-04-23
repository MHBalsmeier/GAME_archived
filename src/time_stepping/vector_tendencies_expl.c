/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, the calculation of the explicit part of the momentum equation is managed.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../game_types.h"
#include "../spatial_operators/spatial_operators.h"
#include "../constituents/constituents.h"
#include "../subgrid_scale/subgrid_scale.h"

int vector_tendencies_expl(State *state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irrev, Config *config, int no_rk_step, double delta_t)
{
	/*
	Managing momentum advection
	---------------------------
	*/
	if (no_rk_step == 1 || config -> totally_first_step_bool == 1)
	{
		// Here, the gaseous flux density is prepared for the generalized Coriolis term.
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			if (config -> assume_lte == 0)
			{
				diagnostics -> scalar_field_placeholder[i] = density_gas(state, i);
			}
			else
			{
				diagnostics -> scalar_field_placeholder[i] = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i];
			}
		}
		scalar_times_vector(diagnostics -> scalar_field_placeholder, state -> wind, diagnostics -> flux_density, grid);
		// Now, the "potential vorticity" is evaluated.
		calc_pot_vort(state -> wind, diagnostics -> scalar_field_placeholder, diagnostics, grid, dualgrid);
		// Now, the generalized Coriolis term is evaluated.
		vorticity_flux(diagnostics -> flux_density, diagnostics -> pot_vort, forcings -> pot_vort_tend, grid, dualgrid);
		// Kinetic energy is prepared for the gradient term of the Lamb transformation.
		inner_product(state -> wind, state -> wind, diagnostics -> v_squared, grid);
		// Taking the gradient of the kinetic energy
		grad(diagnostics -> v_squared, forcings -> v_squared_grad, grid);
    }
    
    /*
    Managing momentum diffusion
    ---------------------------
    */
    if (no_rk_step == 0)
    {
		// momentum diffusion and dissipation (only updated at the first RK step)
		// horizontal momentum diffusion
		if (config -> momentum_diff_h == 1)
		{
			hori_momentum_diffusion(state, diagnostics, irrev, config, grid, dualgrid);
		}
		// vertical momentum diffusion
		if (config -> momentum_diff_v == 1)
		{
			vert_momentum_diffusion(state, diagnostics, irrev, grid, config, delta_t);
		}
		// planetary boundary layer
		if (config -> pbl_scheme > 0)
		{
			pbl_wind_tendency(state, diagnostics, irrev, grid, config, delta_t);
		}
		// calculation of the dissipative heating rate
		if (config -> momentum_diff_h == 1 || config -> pbl_scheme > 0)
		{
			simple_dissipation_rate(state, irrev, grid);
		}
		// Due to condensates, the friction acceleration needs to get a deceleration factor.
		if (config -> assume_lte == 0)
		{
			scalar_times_vector(irrev -> pressure_gradient_decel_factor, irrev -> friction_acc, irrev -> friction_acc, grid);
		}
	}
	
    // Now the explicit forces are added up.
    double old_weight, new_weight;
    new_weight = 1;
    if (no_rk_step == 1)
    {
    	new_weight = 0.5;
    }
	old_weight = 1 - new_weight;
	// the weights for the pressure gradient
	double old_hor_pgrad_weight, current_hor_pgrad_weight, current_ver_pgrad_weight;
	current_hor_pgrad_weight = 0.5 + config -> impl_thermo_weight;
	old_hor_pgrad_weight = 1 - current_hor_pgrad_weight;
	current_ver_pgrad_weight = 1 - config -> impl_thermo_weight;
    int layer_index, h_index;
    #pragma omp parallel for private(layer_index, h_index)
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	layer_index = i/NO_OF_VECTORS_PER_LAYER;
    	h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
    	// upper and lower boundary
        if (i < NO_OF_SCALARS_H || i >= NO_OF_VECTORS - NO_OF_SCALARS_H)
        {
            state_tendency -> wind[i] = 0;
        }
        // horizontal case
        else if (h_index >= NO_OF_SCALARS_H
    	// checking for shading
    	&& NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_vector[h_index - NO_OF_SCALARS_H])
    	{
    		state_tendency -> wind[i] =
    		old_weight*state_tendency -> wind[i] + new_weight*(
    		// explicit component of pressure gradient acceleration
    		// old time step component
    		old_hor_pgrad_weight*forcings -> pgrad_acc_old[i]
    		// current time step component
    		- current_hor_pgrad_weight*(forcings -> pressure_gradient_acc_neg_nl[i] + forcings -> pressure_gradient_acc_neg_l[i])
    		// generalized Coriolis term
    		+ forcings -> pot_vort_tend[i]
    		// kinetic energy term
    		- 0.5*forcings -> v_squared_grad[i]
    		// momentum diffusion
    		+ irrev -> friction_acc[i]);
    	}
        // vertical case
    	else if (h_index < NO_OF_SCALARS_H
    	// checking for shading
    	&& NO_OF_LAYERS - layer_index > grid -> no_of_shaded_points_scalar[h_index])
		{
    		state_tendency -> wind[i] =
    		old_weight*state_tendency -> wind[i] + new_weight*(
    		// explicit component of pressure gradient acceleration
    		// current time step component
    		-current_ver_pgrad_weight*(forcings -> pressure_gradient_acc_neg_nl[i] + forcings -> pressure_gradient_acc_neg_l[i])
    		// generalized Coriolis term
    		+ forcings -> pot_vort_tend[i]
    		// kinetic energy term
    		- 0.5*forcings -> v_squared_grad[i]
    		// momentum diffusion
    		+ irrev -> friction_acc[i]);
		}
    }
    return 0;
}
    
    
    
    
    	
    
    
    
    
    
    
