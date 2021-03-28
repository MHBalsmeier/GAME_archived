/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
*/

#include "../enum_and_typedefs.h"
#include "../settings.h"
#include "../spatial_operators/spatial_operators.h"
#include "time_stepping.h"
#include "../physics/physics.h"
#include "../diagnostics/diagnostics.h"
#include <geos95.h>
#include <stdlib.h>
#include <stdio.h>

int temperature_step(State *, State *, State *, Diagnostics *, Config_info *, double, int);

int manage_rkhevi(State *state_old, State *state_new, Extrapolation_info *extrapolation_info, Grid *grid, Dualgrid *dualgrid, Scalar_field radiation_tendency, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irreversible_quantities, Config_info *config_info, double delta_t, double time_coordinate, int total_step_counter)
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
		// Update of the pressure gradient.
		if (i == 0)
		{
			manage_pressure_gradient(state_new, grid, dualgrid, diagnostics, forcings, extrapolation_info, irreversible_quantities, config_info);
		}
		// Only the horizontal momentum is a forward tendency.
		vector_tendencies_expl(state_new, state_tendency, grid, dualgrid, diagnostics, forcings, irreversible_quantities, config_info, slow_update_bool, i, delta_t);
	    // time stepping for the horizontal momentum can be directly executed
	    
	    for (int j = 0; j < NO_OF_VECTORS; ++j)
	    {
	    	layer_index = j/NO_OF_VECTORS_PER_LAYER;
	    	h_index = j - layer_index*NO_OF_VECTORS_PER_LAYER;
	    	if (h_index >= NO_OF_SCALARS_H)
	    	{
	    		state_new -> velocity_gas[j] = state_old -> velocity_gas[j] + delta_t*state_tendency -> velocity_gas[j];
	    	}
	    }
		// Horizontal velocity can be considered to be updated from now on.

		// 2.) Explicit component of the generalized density equations.
		// ------------------------------------------------------------
	    // Radiation is updated here.
		if (config_info -> rad_on == 1 && config_info -> rad_update == 1 && i == 0)
		{
			printf("Starting update of radiative fluxes ...\n");
			// Fortran needs pointers, this is why this is necessary
			int no_of_scalars = NO_OF_SCALARS;
			int no_of_vectors = NO_OF_VECTORS;
			int no_of_vectors_per_layer = NO_OF_VECTORS_PER_LAYER;
			int no_of_constituents = NO_OF_CONSTITUENTS;
			int no_of_condensed_constituents = NO_OF_CONDENSED_CONSTITUENTS;
			int no_of_layers = NO_OF_LAYERS;
			calc_radiative_flux_convergence(grid -> latitude_scalar, grid -> longitude_scalar, grid -> z_scalar, grid -> z_vector,
			state_new -> mass_densities, state_new -> temperature_gas, radiation_tendency, &no_of_scalars, &no_of_vectors, &no_of_vectors_per_layer, &no_of_layers, &no_of_constituents, 
			&no_of_condensed_constituents, &time_coordinate);
			printf("Update of radiative fluxes completed.\n");
		}
		// Temperature diffusion gets updated here, but only at the first RK step and if heat conduction is switched on.
		if (i == 0 && (config_info -> temperature_diff_h == 1 || config_info -> temperature_diff_v == 1))
		{
			// Now we need to calculate the temperature diffusion coefficients.
		    calc_temp_diffusion_coeffs(state_new, config_info, irreversible_quantities, diagnostics, delta_t, grid);
		    // The diffusion of the temperature depends on its gradient.
			grad(state_new -> temperature_gas, diagnostics -> temperature_gradient, grid);
			// Now the diffusive temperature flux density can be obtained.
		    scalar_times_vector_scalar_h_v(irreversible_quantities -> scalar_diffusion_coeff_numerical_h, irreversible_quantities -> scalar_diffusion_coeff_numerical_v,
		    diagnostics -> temperature_gradient, diagnostics -> flux_density, grid);
		    // The divergence of the diffusive temperature flux density is the diffusive temperature heating.
		    divv_h(diagnostics -> flux_density, irreversible_quantities -> temperature_diffusion_heating, grid);
		    add_vertical_divv(diagnostics -> flux_density, irreversible_quantities -> temperature_diffusion_heating, grid);
		}
		if (i == 0)
		{
			scalar_tendencies_expl(state_new, state_tendency, grid, dualgrid, delta_t, radiation_tendency, diagnostics, forcings, irreversible_quantities, config_info, i, state_old -> velocity_gas);
		}
		if (i == 1)
		{	
			scalar_tendencies_expl(state_new, state_tendency, grid, dualgrid, delta_t, radiation_tendency, diagnostics, forcings, irreversible_quantities, config_info, i, state_new -> velocity_gas);
		}
		
		// 3.) A pre-conditioned new temperature field, only containing explicit entropy and mass density tendencies (including diabatic forcings).
		// ----------------------------------------------------------------------------------------------------------------------------------------
		temperature_step(state_old, state_new, state_tendency, diagnostics, config_info, delta_t, 0);

		// 4.) Vertical sound wave solver.
		// -------------------------------
		three_band_solver_ver_waves(state_old, state_new, state_tendency, diagnostics, config_info, delta_t, grid, i);
		// Vertical velocity can be seen as updated from now on.
		
		// 5.) Solving the implicit component of the generalized density equations for tracers.
		// ------------------------------------------------------------------------------------
		if (NO_OF_CONSTITUENTS > 1)
		{
			three_band_solver_gen_densitites(state_old, state_new, state_tendency, diagnostics, config_info, delta_t, grid);
		}
    }

	// in this case, a large time step has been taken, which we modify into a small step here    
    if (slow_update_bool == 1 && config_info -> adv_sound_ratio > 1)
    {
    	linear_combine_two_states(state_old, state_new, state_new, 1 - delta_t_small/delta_t, delta_t_small/delta_t);
    	// this is for thermodynamic consistency
		temperature_step(state_old, state_new, state_tendency, diagnostics, config_info, delta_t_small, 1);
    }
    
    return 0;
}

int temperature_step(State *state_old, State *state_new, State *state_tendency, Diagnostics *diagnostics, Config_info *config_info, double delta_t, int write2new)
{
	// explicit temperature diagnostics based on linearization of internal energy
    double nominator, denominator, entropy_density_gas_0, entropy_density_gas_1, density_gas_0, density_gas_1, delta_density_gas, delta_entropy_density, temperature_0, specific_entropy_gas_0, specific_entropy_gas_1, c_g_v, c_g_p;
    double beta = get_impl_thermo_weight();
    double alpha = 1 - beta;
	#pragma omp parallel for private(nominator, denominator, entropy_density_gas_0, entropy_density_gas_1, density_gas_0, density_gas_1, delta_density_gas, delta_entropy_density, temperature_0, specific_entropy_gas_0, specific_entropy_gas_1, c_g_v, c_g_p)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	// Difference of the mass densities of the gas phase.
    	density_gas_0 = 0;
    	density_gas_1 = 0;
		int no_of_relevant_constituents = 0;
		if (config_info -> assume_lte == 0)
		{
			no_of_relevant_constituents = NO_OF_GASEOUS_CONSTITUENTS;
		}
		if (config_info -> assume_lte == 1)
		{
			no_of_relevant_constituents = 1;
		}
    	for (int j = 0; j < no_of_relevant_constituents; ++j)
    	{
			density_gas_0 += state_old -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + j)*NO_OF_SCALARS + i];
			if (write2new == 0)
			{
				density_gas_1 += state_old -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + j)*NO_OF_SCALARS + i]
				+ delta_t*state_tendency -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + j)*NO_OF_SCALARS + i];
			}
			if (write2new == 1)
			{
				density_gas_1 += state_new -> mass_densities[(NO_OF_CONDENSED_CONSTITUENTS + j)*NO_OF_SCALARS + i];
			}
    	}
    	delta_density_gas = density_gas_1 - density_gas_0;
    	
    	entropy_density_gas_0 = 0;
    	entropy_density_gas_1 = 0;
    	for (int j = 0; j < no_of_relevant_constituents; ++j)
    	{
			entropy_density_gas_0 += state_old -> entropy_densities[j*NO_OF_SCALARS + i];
			if (write2new == 0)
			{
				entropy_density_gas_1 += state_old -> entropy_densities[j*NO_OF_SCALARS + i]
				+ delta_t*state_tendency -> entropy_densities[j*NO_OF_SCALARS + i];
			}
			if (write2new == 1)
			{
				entropy_density_gas_1 += state_new -> entropy_densities[j*NO_OF_SCALARS + i];
			}
    	}
    	delta_entropy_density = entropy_density_gas_1 - entropy_density_gas_0;
    	
    	// Specific entropies of the gas phase of the two time steps.
    	specific_entropy_gas_0 = entropy_density_gas_0/density_gas_0;
    	specific_entropy_gas_1 = entropy_density_gas_1/density_gas_1;
    	
    	// The temperature of the gas phase of the old time step.
    	temperature_0 = state_old -> temperature_gas[i];
    	
		// Determining the thermodynamic properties of the gas phase.
    	c_g_v = spec_heat_cap_diagnostics_v(state_old, i, config_info);
    	c_g_p = spec_heat_cap_diagnostics_p(state_old, i, config_info);
    	
    	nominator = c_g_v*density_gas_0*temperature_0 + (alpha*c_g_p*temperature_0 - alpha*specific_entropy_gas_0*temperature_0)*delta_density_gas + alpha*temperature_0*delta_entropy_density;
    	denominator = c_g_v*density_gas_0 + (c_g_v + beta*specific_entropy_gas_1 - beta*c_g_p)*delta_density_gas - beta*delta_entropy_density;
    	if (write2new == 0)
    	{
			diagnostics -> temperature_gas_explicit[i] = nominator/denominator;
		}
    	if (write2new == 1)
    	{
			state_new -> temperature_gas[i] = nominator/denominator;
		}
    }
    return 0;
}






