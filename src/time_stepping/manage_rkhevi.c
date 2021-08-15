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

int create_rad_array_scalar(double [], double [], int);
int create_rad_array_scalar_h(double [], double [], int);
int create_rad_array_mass_den(double [], double [], int);
int create_rad_array_vector(double [], double [], int);
int remap_to_original(double [], double [], int);
int remap_to_original_scalar_h(double [], double [], int);

int manage_rkhevi(State *state_old, State *state_new, Soil *soil, Grid *grid, Dualgrid *dualgrid, Radiation *radiation, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irrev, Config_info *config_info, double delta_t, double time_coordinate, int total_step_counter)
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
       
    // interaction with soil (only useful if real radiation is on)
    if (config_info -> rad_on == 1)
    {
    	soil_interaction(soil, diagnostics, radiation, delta_t);
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
		
		// 0.) diagnosing the temperature
		temperature_diagnostics(state_new, grid, diagnostics);
		
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
	    
	    #pragma omp parallel for private(layer_index, h_index)
	    for (int j = 0; j < NO_OF_VECTORS; ++j)
	    {
	    	layer_index = j/NO_OF_VECTORS_PER_LAYER;
	    	h_index = j - layer_index*NO_OF_VECTORS_PER_LAYER;
	    	if (h_index >= NO_OF_SCALARS_H)
	    	{
	    		state_new -> wind[j] = state_old -> wind[j] + delta_t*state_tendency -> wind[j];
	    	}
	    }
		// Horizontal velocity can be considered to be updated from now on.

		// 2.) Explicit component of the generalized density equations.
		// ------------------------------------------------------------
	    // Radiation is updated here.
		if (config_info -> rad_on > 0 && config_info -> rad_update == 1 && i == 0)
		{
			printf("Starting update of radiative fluxes ...\n");
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
				create_rad_array_scalar_h(soil -> temperature, radiation -> temp_sfc_rad, rad_block_index);
				create_rad_array_scalar(grid -> z_scalar, radiation -> z_scal_rad, rad_block_index);
				create_rad_array_vector(grid -> z_vector, radiation -> z_vect_rad, rad_block_index);
				create_rad_array_mass_den(state_old -> rho, radiation -> rho_rad, rad_block_index);
				create_rad_array_scalar(diagnostics -> temperature_gas, radiation -> temp_rad, rad_block_index);
				// calling the radiation routine
				// RTE+RRTMGP
				if (config_info -> rad_on == 1)
				{
					calc_radiative_flux_convergence(radiation -> lat_scal_rad,
					radiation -> lon_scal_rad,
					radiation -> z_scal_rad,
					radiation -> z_vect_rad,
					radiation -> rho_rad,
					radiation -> temp_rad,
					radiation -> rad_tend_rad,
					radiation -> temp_sfc_rad,
					radiation -> sfc_sw_in_rad,
					radiation -> sfc_lw_out_rad,
					&no_of_scalars, &no_of_layers,
					&no_of_constituents, &no_of_condensed_constituents,
					&time_coordinate);
				}
				// Held-Suarez
				if (config_info -> rad_on == 2)
				{
					held_suar(radiation -> lat_scal_rad, radiation -> z_scal_rad, radiation -> rho_rad, radiation -> temp_rad, radiation -> rad_tend_rad);
				}
				// filling the actual radiation tendency
				remap_to_original(radiation -> rad_tend_rad, radiation -> radiation_tendency, rad_block_index);
				remap_to_original_scalar_h(radiation -> sfc_sw_in_rad, radiation -> sfc_sw_in, rad_block_index);
				remap_to_original_scalar_h(radiation -> sfc_lw_out_rad, radiation -> sfc_lw_out, rad_block_index);
			}
			printf("Update of radiative fluxes completed.\n");
		}
		scalar_tendencies_expl(state_new, state_tendency, soil, grid, dualgrid, delta_t, radiation -> radiation_tendency, diagnostics, forcings, irrev, config_info, i);

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

	// in this case, a large time step has been taken, which we modify into a small step here
    if (slow_update_bool == 1 && config_info -> adv_sound_ratio > 1)
    {
    	linear_combine_two_states(state_old, state_new, state_new, 1 - delta_t_small/delta_t, delta_t_small/delta_t, grid);
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


int remap_to_original_scalar_h(double in[], double out[], int rad_block_index)
{
	/*
	reverses what create_rad_array_scalar_h has done
	*/
	// loop over all elements of the resulting array
	for (int i = 0; i < NO_OF_SCALARS_RAD_PER_LAYER; ++i)
	{
		out[rad_block_index*NO_OF_SCALARS_RAD_PER_LAYER + i] = in[i];
	}
	return 0;
}






