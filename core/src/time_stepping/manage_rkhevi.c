/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
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
#include <geos95.h>
#include <stdlib.h>
#include <stdio.h>

int temperature_step(State *, State *, State *, Diagnostics *, Config_info *, double, int);
int create_rad_array_scalar(double [], double [], int);
int create_rad_array_scalar_h(double [], double [], int);
int create_rad_array_mass_den(double [], double [], int);
int create_rad_array_vector(double [], double [], int);
int remap_to_original(double [], double [], int);

int manage_rkhevi(State *state_old, State *state_new, Extrapolation_info *extrapolation_info, Grid *grid, Dualgrid *dualgrid, Radiation *radiation, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irrev, Config_info *config_info, double delta_t, double time_coordinate, int total_step_counter)
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
			manage_pressure_gradient(state_new, grid, dualgrid, diagnostics, forcings, extrapolation_info, irrev, config_info);
		}
		// Only the horizontal momentum is a forward tendency.
		vector_tendencies_expl(state_new, state_tendency, grid, dualgrid, diagnostics, forcings, irrev, config_info, slow_update_bool, i, delta_t);
	    // time stepping for the horizontal momentum can be directly executed
	    
	    #pragma omp parallel for private(layer_index, h_index)
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
		if (config_info -> rad_on > 0 && config_info -> rad_update == 1 && i == 0)
		{
			printf("Starting update of radiative fluxes ...\n");
			if (config_info -> rad_on == 1)
			{
				int no_of_scalars = NO_OF_SCALARS_RAD;
				int no_of_constituents = NO_OF_CONSTITUENTS;
				int no_of_condensed_constituents = NO_OF_CONDENSED_CONSTITUENTS;
				int no_of_layers = NO_OF_LAYERS;
				// loop over all radiation blocks
				for (int rad_block_index = 0; rad_block_index < NO_OF_RAD_BLOCKS; ++rad_block_index)
				{
					// remapping all the arrays
					create_rad_array_scalar_h(grid -> latitude_scalar, radiation -> lat_scal_rad, rad_block_index);
					create_rad_array_scalar_h(grid -> longitude_scalar, radiation -> lon_scal_rad, rad_block_index);
					create_rad_array_scalar(grid -> z_scalar, radiation -> z_scal_rad, rad_block_index);
					create_rad_array_vector(grid -> z_vector, radiation -> z_vect_rad, rad_block_index);
					create_rad_array_mass_den(state_old -> mass_densities, radiation -> mass_den_rad, rad_block_index);
					create_rad_array_scalar(state_old -> temperature_gas, radiation -> temp_rad, rad_block_index);
					// calling the radiation routine
					calc_radiative_flux_convergence(radiation -> lat_scal_rad,
					radiation -> lon_scal_rad,
					radiation -> z_scal_rad,
					radiation -> z_vect_rad,
					radiation -> mass_den_rad,
					radiation -> temp_rad,
					radiation -> rad_tend_rad,
					&no_of_scalars, &no_of_layers,
					&no_of_constituents, &no_of_condensed_constituents,
					&time_coordinate);
					// filling the actual radiation tendency
					remap_to_original(radiation -> rad_tend_rad, radiation -> radiation_tendency, rad_block_index);
				}
			}
			if (config_info -> rad_on == 2)
			{
				held_suar(grid -> latitude_scalar, grid -> z_scalar, state_old -> mass_densities, state_old -> temperature_gas, radiation -> radiation_tendency);
			}
			printf("Update of radiative fluxes completed.\n");
		}
		if (i == 0)
		{
			scalar_tendencies_expl(state_new, state_tendency, grid, dualgrid, delta_t, radiation -> radiation_tendency, diagnostics, forcings, irrev, config_info, i, state_old -> velocity_gas);
		}
		if (i == 1)
		{	
			scalar_tendencies_expl(state_new, state_tendency, grid, dualgrid, delta_t, radiation -> radiation_tendency, diagnostics, forcings, irrev, config_info, i, state_new -> velocity_gas);
		}
		
		// 3.) A pre-conditioned new temperature field, only containing explicit entropy and mass density tendencies (including diabatic forcings).
		// ----------------------------------------------------------------------------------------------------------------------------------------
		temperature_step(state_old, state_new, state_tendency, diagnostics, config_info, delta_t, 0);

		// 4.) Vertical sound wave solver.
		// -------------------------------
		three_band_solver_ver_waves(state_old, state_new, state_tendency, diagnostics, config_info, delta_t, grid, i);
		// Vertical velocity can be seen as updated from now on.
		// this is for stability
		temperature_step(state_old, state_new, state_tendency, diagnostics, config_info, delta_t, 1);
		
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
		// this is for stability
		temperature_step(state_old, state_new, state_tendency, diagnostics, config_info, delta_t_small, 1);
    }
	
	// nesting
	if (config_info -> regional_switch == 1)
	{
		bc_setter();
	}
    
    return 0;
}

int temperature_step(State *state_old, State *state_new, State *state_tendency, Diagnostics *diagnostics, Config_info *config_info, double delta_t, int write2new)
{
	// temperature step based on linearization of internal energy
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
    	
		// determining the thermodynamic properties of the gas phase
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


int create_rad_array_scalar(double in[], double out[], int rad_block_index)
{
	/*
	cuts out a slice of a scalar field for hand-over to the radiation routine (done for RAM efficiency reasons)
	*/
	int layer_index, h_index;
	// loop over all elements of the resulting array
	for (int i = 0; i < NO_OF_SCALARS_RAD; ++i)
	{
		layer_index = i/NO_OF_SCALARS_RAD_PER_LAYER;
		h_index = i - layer_index*NO_OF_SCALARS_RAD_PER_LAYER;
		out[i] = in[rad_block_index*NO_OF_SCALARS_RAD_PER_LAYER + h_index + layer_index*NO_OF_SCALARS_H];
	}
	return 0;
}

int create_rad_array_scalar_h(double in[], double out[], int rad_block_index)
{
	/*
	cuts out a slice of a horizontal scalar field for hand-over to the radiation routine (done for RAM efficiency reasons)
	*/
	// loop over all elements of the resulting array
	for (int i = 0; i < NO_OF_SCALARS_RAD_PER_LAYER; ++i)
	{
		out[i] = in[rad_block_index*NO_OF_SCALARS_RAD_PER_LAYER + i];
	}
	return 0;
}

int create_rad_array_mass_den(double in[], double out[], int rad_block_index)
{
	/*
	same thing as create_rad_array_scalar, only for a mass density field
	*/
	int layer_index, h_index;
	for (int const_id = 0; const_id < NO_OF_CONSTITUENTS; ++const_id)
	{
		// loop over all elements of the resulting array
		for (int i = 0; i < NO_OF_SCALARS_RAD; ++i)
		{
			layer_index = i/NO_OF_SCALARS_RAD_PER_LAYER;
			h_index = i - layer_index*NO_OF_SCALARS_RAD_PER_LAYER;
			out[const_id*NO_OF_SCALARS_RAD + i]
			= in[const_id*NO_OF_SCALARS + rad_block_index*NO_OF_SCALARS_RAD_PER_LAYER + h_index + layer_index*NO_OF_SCALARS_H];
		}
	}
	return 0;
}

int create_rad_array_vector(double in[], double out[], int rad_block_index)
{
	/*
	cuts out a slice of a vector field for hand-over to the radiation routine (done for RAM efficiency reasons),
	only the vertical vector points are taken into account since only they are needed by the radiation
	*/
	int layer_index, h_index;
	// loop over all elements of the resulting array
	for (int i = 0; i < NO_OF_SCALARS_RAD + NO_OF_SCALARS_RAD_PER_LAYER; ++i)
	{
		layer_index = i/NO_OF_SCALARS_RAD_PER_LAYER;
		h_index = i - layer_index*NO_OF_SCALARS_RAD_PER_LAYER;
		out[i] = in[rad_block_index*NO_OF_SCALARS_RAD_PER_LAYER + h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
	}
	return 0;
}

int remap_to_original(double in[], double out[], int rad_block_index)
{
	/*
	reverses what create_rad_array_scalar has done
	*/
	int layer_index, h_index;
	// loop over all elements of the resulting array
	for (int i = 0; i < NO_OF_SCALARS_RAD; ++i)
	{
		layer_index = i/NO_OF_SCALARS_RAD_PER_LAYER;
		h_index = i - layer_index*NO_OF_SCALARS_RAD_PER_LAYER;
		out[rad_block_index*NO_OF_SCALARS_RAD_PER_LAYER + h_index + layer_index*NO_OF_SCALARS_H] = in[i];
	}
	return 0;
}



