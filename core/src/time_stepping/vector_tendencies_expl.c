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
#include "../thermodynamics/thermodynamics.h"

int vector_tendencies_expl(State *state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irrev, Config_info *config_info, int update_advection, int no_rk_step, double delta_t)
{
	// momentum advection
	if ((update_advection == 1 && no_rk_step == 1) || config_info -> totally_first_step_bool == 1)
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
    
    // momentum diffusion and dissipation (only updated at the first RK step and if advection is updated as well)
    if (no_rk_step == 0 && update_advection == 1)
    {
    	// horizontal momentum diffusion
    	if (config_info -> momentum_diff_h == 1)
    	{
			hori_momentum_diffusion(state, diagnostics, irrev, config_info, grid, dualgrid, delta_t);
		}
		// vertical momentum diffusion
		if (config_info -> momentum_diff_v == 1)
		{
			vert_momentum_diffusion(state, diagnostics, irrev, grid, config_info, delta_t);
		}
		// This is an explicit friction ansatz in the boundary layer, comparable to what is required in the Held-Suarez test.
		if (config_info -> explicit_boundary_layer == 1)
		{
			// some parameters
			double bndr_lr_height = 1e3; // boundary layer height
			double bndr_lr_visc_max = 1.0/86400; // maximum friction coefficient in the boundary layer
			double e_folding_height = 0.5*bndr_lr_height;
			double z_agl;
			int layer_index, h_index, vector_index;
			#pragma omp parallel for private(layer_index, h_index, vector_index, z_agl)
			for (int i = 0; i < NO_OF_H_VECTORS; ++i)
			{
				layer_index = i/NO_OF_VECTORS_H;
				h_index = i - layer_index*NO_OF_VECTORS_H;
				vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
				// height above ground level
				z_agl = grid -> z_vector[vector_index]
				- 0.5*(grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> from_index[h_index]]
				+ grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> to_index[h_index]]);
				// adding the boundary layer friction
				if (z_agl < bndr_lr_height)
				{
					irrev -> friction_acc[vector_index]
					+= -bndr_lr_visc_max*(exp(-z_agl/e_folding_height) - exp(-bndr_lr_height/e_folding_height))
					/(1 - exp(-bndr_lr_height/e_folding_height))
					*state -> velocity_gas[vector_index];
				}
			}
		}
		// calculation of the dissipative heating rate
		if (config_info -> dissipative_heating == 1)
		{
			simple_dissipation_rate(state, irrev, grid);
		}
		// Due to condensates, the friction acceleration needs to get a deceleration factor.
		if (config_info -> assume_lte == 0)
		{
			scalar_times_vector(irrev -> pressure_gradient_decel_factor, irrev -> friction_acc, irrev -> friction_acc, grid);
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
    		+ irrev -> friction_acc[i]);
		}
    }
    return 0;
}
    
    
    
    
    	
    
    
    
    
    
    
