/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/
/*
In this source file, the calculation of the explicit part of the momentum equation is managed.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../../enum_and_typedefs.h"
#include "../../settings.h"
#include "../../spatial_operators/spatial_operators.h"
#include "atmostracers.h"
#include "../../diagnostics/diagnostics.h"

int integrate_momentum(State *state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irreversible_quantities, Config_info *config_info)
{
	// Here, the gaseous flux density is prepared for the generalized Coriolis term.
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = density_gas(state, i);
	}
    scalar_times_vector(diagnostics -> scalar_field_placeholder, state -> velocity_gas, diagnostics -> flux_density, grid);
    // Now, the potential vorticity is evaluated.
    calc_pot_vort(state -> velocity_gas, diagnostics -> scalar_field_placeholder, diagnostics, grid, dualgrid);
    // Now, the generalized Coriolis term is evaluated.
    vorticity_flux(diagnostics -> flux_density, diagnostics -> pot_vort, forcings -> pot_vort_tend, grid, dualgrid);
    // Horizontal kinetic energy is prepared for the gradient term of the Lamb transformation.
    kinetic_energy(state -> velocity_gas, diagnostics -> e_kin_h, grid, 0);
    // taking the gradient of the horizontal kinetic energy
    grad(diagnostics -> e_kin_h, forcings -> e_kin_h_grad, grid);
    // Now the explicit forces are added up.
    int layer_index, h_index;
    double metric_term, vertical_velocity, hor_non_trad_cori_term, expl_pgrad_weight;
	expl_pgrad_weight = 1 - get_impl_thermo_weight();
    #pragma omp parallel for private(layer_index, h_index, metric_term, vertical_velocity, hor_non_trad_cori_term)
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	layer_index = i/NO_OF_VECTORS_PER_LAYER;
    	h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
    	// upper and lower boundary
        if (i < NO_OF_SCALARS_H || i >= NO_OF_VECTORS - NO_OF_SCALARS_H)
        {
            state_tendency -> velocity_gas[i] = 0;
        }
        else
        {
        	// horizontal case
    		if (h_index >= NO_OF_SCALARS_H
        	// checking for shading
        	&& NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_vector[h_index - NO_OF_SCALARS_H])
    		{
    			// determining w at the edge
    			remap_verpri2horpri_vector(state -> velocity_gas, layer_index, h_index - NO_OF_SCALARS_H, &vertical_velocity, grid);
    			// deep atmospehre metric term -w/r*v_h
    			metric_term = -vertical_velocity*state -> velocity_gas[i]/(RADIUS + grid -> z_vector[i]);
    			// (-f_y*w*i)*n = -f_y*w*(i*n) = -f_y*w*cos(direction)
    			hor_non_trad_cori_term = -vertical_velocity*dualgrid -> f_vec[2*NO_OF_VECTORS_H + h_index - NO_OF_SCALARS_H];
        		state_tendency -> velocity_gas[i] =
        		// full pressure gradient acceleration
        		forcings -> pressure_gradient_acc[i]
        		// generalized Coriolis term
        		+ forcings -> pot_vort_tend[i]
        		// gravity
        		- grid -> gravity_m[i]
        		// horizontal kinetic energy term
        		- forcings -> e_kin_h_grad[i]
        		// the term that results from the horizontal Coriolis vector
        		+ hor_non_trad_cori_term
        		// deep atmosphere metric term
        		+ metric_term
        		// momentum diffusion
        		+ irreversible_quantities -> friction_acc[i];
    		}
    		// vertical case
        	if (h_index < NO_OF_SCALARS_H
        	// checking for shading
        	&& NO_OF_LAYERS - layer_index > grid -> no_of_shaded_points_scalar[h_index])
        	{
        		state_tendency -> velocity_gas[i] =
        		// generalized Coriolis term
        		forcings -> pot_vort_tend[i]
        		// horizontal kinetic energy term
        		- forcings -> e_kin_h_grad[i]
        		// gravity
        		- grid -> gravity_m[i]
        		// explicit part of the pressure gradient acceleration
        		+ expl_pgrad_weight*forcings -> pressure_gradient_acc[i]
        		// momentum diffusion
        		+ irreversible_quantities -> friction_acc[i];
    		}
        }
    }
    return 0;
}
    
    
    
    
    	
    
    
    
    
    
    
