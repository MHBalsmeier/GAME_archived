/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/
/*
In this source file, the calculation of the explicit part of the momentum equation is managed.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../enum_and_typedefs.h"
#include "../settings.h"
#include "../spatial_operators/spatial_operators.h"
#include "../diagnostics/diagnostics.h"

int vector_tendencies_expl(State *state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irreversible_quantities, Config_info *config_info, int update_advection, int no_rk_step, double delta_t)
{
	// momentum advection
	if (update_advection == 1)
	{
		// Here, the gaseous flux density is prepared for the generalized Coriolis term.
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			if (config_info -> assume_lte == 0)
			{
				diagnostics -> scalar_field_placeholder[i] = density_gas(state, i);
			}
			else
			{
				diagnostics -> scalar_field_placeholder[i] = state -> mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i];
			}
		}
		scalar_times_vector(diagnostics -> scalar_field_placeholder, state -> velocity_gas, diagnostics -> flux_density, grid);
		// Now, the potential vorticity is evaluated.
		calc_pot_vort(state -> velocity_gas, diagnostics -> scalar_field_placeholder, diagnostics, grid, dualgrid);
		// Now, the generalized Coriolis term is evaluated.
		vorticity_flux(diagnostics -> flux_density, diagnostics -> pot_vort, forcings -> pot_vort_tend, grid, dualgrid);
		// Kinetic energy is prepared for the gradient term of the Lamb transformation.
		kinetic_energy(state -> velocity_gas, diagnostics -> e_kin, grid);
		// Taking the gradient of the kinetic energy
		grad(diagnostics -> e_kin, forcings -> e_kin_grad, grid);
    }
    // momentum diffusion and dissipation
    if (no_rk_step == 0 && config_info -> momentum_diff_h == 1 && update_advection == 1)
    {
		momentum_diff_diss(state, diagnostics, irreversible_quantities, config_info, grid, dualgrid, delta_t);
		// Due to condensates, the friction acceleration needs to get a deceleration factor.
		if (config_info -> assume_lte == 0)
		{
			scalar_times_vector(irreversible_quantities -> pressure_gradient_decel_factor, irreversible_quantities -> friction_acc, irreversible_quantities -> friction_acc, grid);
		}
    }
    // Now the explicit forces are added up.
    int layer_index, h_index;
    double old_weight, new_weight;
    new_weight = 1;
    if (no_rk_step == 1)
    {
    	new_weight = 0.5;
    }
	old_weight = 1 - new_weight;
    #pragma omp parallel for private(layer_index, h_index)
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	layer_index = i/NO_OF_VECTORS_PER_LAYER;
    	h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
    	// upper and lower boundary
        if (i < NO_OF_SCALARS_H || i >= NO_OF_VECTORS - NO_OF_SCALARS_H)
        {
            state_tendency -> velocity_gas[i] = 0;
        }
        else if ((h_index >= NO_OF_SCALARS_H
    	// checking for shading
    	&& NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_vector[h_index - NO_OF_SCALARS_H])
    	|| (h_index < NO_OF_SCALARS_H
    	// checking for shading
    	&& NO_OF_LAYERS - layer_index > grid -> no_of_shaded_points_scalar[h_index]))
		{
    		state_tendency -> velocity_gas[i] =
    		old_weight*state_tendency -> velocity_gas[i] + new_weight*(
    		// explicit component of pressure gradient acceleration
    		+ forcings -> pressure_gradient_acc_expl[i]
    		// generalized Coriolis term
    		+ forcings -> pot_vort_tend[i]
    		// kinetic energy term
    		- forcings -> e_kin_grad[i]
    		// gravity
    		- grid -> gravity_m[i]
    		// momentum diffusion
    		+ irreversible_quantities -> friction_acc[i]);
		}
    }
    return 0;
}
    
    
    
    
    	
    
    
    
    
    
    
