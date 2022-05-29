/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This is the horizontal (explicit) part of the constituent integration.
*/

#include <stdlib.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../spatial_operators/spatial_operators.h"
#include "../constituents/constituents.h"
#include "../subgrid_scale/subgrid_scale.h"

int scalar_tendencies_expl(State *state_old, State *state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Diagnostics *diagnostics, Forcings *forcings,
Irreversible_quantities *irrev, Config *config, int rk_step)
{
	/*
	This function manages the calculation of the explicit part of the scalar tendencies.
	*/
	
	/*
	Firstly, some things need to prepared.
	--------------------------------------
	*/
	// declaring needed variables
    int scalar_shift_index, scalar_index;
    
    // determining the RK weights
    double old_weight[NO_OF_CONSTITUENTS];
    double new_weight[NO_OF_CONSTITUENTS];
    for (int i = 0; i < NO_OF_CONSTITUENTS; ++i)
    {
		new_weight[i] = 1.0;
		if (rk_step == 1 && i != NO_OF_CONDENSED_CONSTITUENTS)
		{
			new_weight[i] = 0.5;
		}
		old_weight[i] = 1.0 - new_weight[i];
    }
    
    // updating the scalar diffusion coefficient if required
    if (rk_step == 0 && (config -> mass_diff_h == 1 || config -> temperature_diff_h == 1))
    {
	    scalar_diffusion_coeffs(state, config, irrev, diagnostics, delta_t, grid, dualgrid);
    }
    
	// Temperature diffusion gets updated at the first RK step if required.
	if (config -> temperature_diff_h == 1 && rk_step == 0)
	{
	    // The diffusion of the temperature depends on its gradient.
		grad(diagnostics -> temperature_gas, diagnostics -> vector_field_placeholder, grid);
		// Now the diffusive temperature flux density can be obtained.
	    scalar_times_vector_h(irrev -> temp_diffusion_coeff_numerical_h, diagnostics -> vector_field_placeholder, diagnostics -> flux_density, grid);
	    // The divergence of the diffusive temperature flux density is the diffusive temperature heating.
	    divv_h(diagnostics -> flux_density, irrev -> temperature_diffusion_heating, grid);
    	// vertical temperature diffusion
	    if (config -> temperature_diff_v == 1)
	    {
	    	scalar_times_vector_v(irrev -> temp_diffusion_coeff_numerical_v, diagnostics -> vector_field_placeholder, diagnostics -> flux_density, grid);
	    	add_vertical_divv(diagnostics -> flux_density, irrev -> temperature_diffusion_heating, grid);
		}
	}
	
	/*
	Now, the actual scalar tendencies can be computed.
	--------------------------------------------------
	*/
	// loop over all constituents
	for (int i = 0; i < NO_OF_CONSTITUENTS; ++i)
	{
		scalar_shift_index = i*NO_OF_SCALARS;
		
        // This is the mass advection, which needs to be carried out for all constituents.
        // -------------------------------------------------------------------------------
        // moist air
		if (i == NO_OF_CONDENSED_CONSTITUENTS)
		{
			scalar_times_vector_h(&state -> rho[scalar_shift_index], state -> wind, diagnostics -> flux_density, grid);
    		divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
		}
		// all other constituents
		else
		{
			scalar_times_vector_h_upstream(&state -> rho[scalar_shift_index], state -> wind, diagnostics -> flux_density, grid);
    		divv_h_tracer(diagnostics -> flux_density, &state -> rho[scalar_shift_index], state -> wind, diagnostics -> flux_density_divv, grid);
		}
		
		// mass diffusion
		if (config -> mass_diff_h == 1)
		{
    		// The diffusion of the tracer density depends on its gradient.
			grad(&state -> rho[scalar_shift_index], diagnostics -> vector_field_placeholder, grid);
			// Now the diffusive mass flux density can be obtained.
			scalar_times_vector_h(irrev -> mass_diffusion_coeff_numerical_h, diagnostics -> vector_field_placeholder, diagnostics -> vector_field_placeholder, grid);
	    	// The divergence of the diffusive mass flux density is the diffusive mass source rate.
			divv_h(diagnostics -> vector_field_placeholder, diagnostics -> scalar_field_placeholder, grid);
			// vertical mass diffusion
			if (config -> mass_diff_v == 1)
			{
				scalar_times_vector_v(irrev -> mass_diffusion_coeff_numerical_v, diagnostics -> vector_field_placeholder, diagnostics -> vector_field_placeholder, grid);
				add_vertical_divv(diagnostics -> vector_field_placeholder, diagnostics -> scalar_field_placeholder, grid);
			}
		}
		
		// adding the tendencies in all grid boxes
		#pragma omp parallel for private(scalar_index)
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
			scalar_index = scalar_shift_index + j;
			state_tendency -> rho[scalar_index]
			= old_weight[i]*state_tendency -> rho[scalar_index]
			+ new_weight[i]*(
			// the advection
			-diagnostics -> flux_density_divv[j]
			// the diffusion
			+ config -> mass_diff_h*diagnostics -> scalar_field_placeholder[j]);
	    }
	    
	    /*
		Explicit component of the rho*theta_v integration
		-------------------------------------------------
		*/
		if (i == NO_OF_CONDENSED_CONSTITUENTS)
		{
			// determining the virtual potential temperature
			#pragma omp parallel for
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				diagnostics -> scalar_field_placeholder[j] = state -> rhotheta_v[j]/state -> rho[scalar_shift_index + j];
			}
			scalar_times_vector_h(diagnostics -> scalar_field_placeholder, diagnostics -> flux_density, diagnostics -> flux_density, grid);
			divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
			// adding the tendencies in all grid boxes
			#pragma omp parallel for private(scalar_index)
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				scalar_index = scalar_shift_index + j;
				state_tendency -> rhotheta_v[j]
				= old_weight[i]*state_tendency -> rhotheta_v[j]
				+ new_weight[i]*(
				// the advection (resolved transport)
				-diagnostics -> flux_density_divv[j]
				// the diabatic forcings
				// weighting factor
				+ 1.0*(
				// dissipation through molecular + turbulent momentum diffusion
				irrev -> heating_diss[j]
				// molecular + turbulent heat transport
				+ irrev -> temperature_diffusion_heating[j]
				// radiation
				+ forcings -> radiation_tendency[j]
				// phase transitions
				+ irrev -> phase_trans_heating_rate[j]
				// this has to be divided by c_p*exner
				)/(spec_heat_capacities_p_gas(0)*(grid -> exner_bg[j] + state -> exner_pert[j])));
			}
		}
	}
	
	return 0;
}

int moisturizer(State *state, double delta_t, Diagnostics *diagnostics, Irreversible_quantities *irrev, Config *config, Grid *grid)
{
	/*
	This function manages the calculation of the phase transition rates.
	*/
	
	// Only if we have multiple constituents, moisture needs to be included.
	if (NO_OF_CONSTITUENTS > 1)
	{
		// calculating the source rates
	    calc_h2otracers_source_rates(state, diagnostics, grid, config, irrev, 2.0*delta_t);
	    int scalar_shift_index, scalar_index;
	    // loop over all constituents
		for (int i = 0; i < NO_OF_CONSTITUENTS; ++i)
		{
			scalar_shift_index = i*NO_OF_SCALARS;
			// the main gaseous constituent has no source rates
			if (i != NO_OF_CONDENSED_CONSTITUENTS)
			{
				#pragma omp parallel for private(scalar_index)
				for (int j = 0; j < NO_OF_SCALARS; ++j)
				{
					scalar_index = scalar_shift_index + j;
					if (i < NO_OF_CONDENSED_CONSTITUENTS)
					{
						state -> rho[scalar_index] = state -> rho[scalar_index] + delta_t*irrev -> mass_source_rates[scalar_index];
					}
					// for the gaseous constituents (apart from the main one), an index shift is necessary
					else
					{
						state -> rho[scalar_index] = state -> rho[scalar_index] + delta_t*irrev -> mass_source_rates[(i - 1)*NO_OF_SCALARS + j];
					}
				}
			}
		}
	}
	return 0;	
}






