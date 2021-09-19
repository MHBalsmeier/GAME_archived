/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This is the horizontal (explicit) part of the constituent integration.
*/

#include "../enum_and_typedefs.h"
#include "atmostracers.h"
#include "../settings.h"
#include "../spatial_operators/spatial_operators.h"
#include "../thermodynamics/thermodynamics.h"
#include "../subgrid_scale/subgrid_scale.h"
#include "stdio.h"
#include "stdlib.h"

int scalar_tendencies_expl(State *state_old, State *state, State *state_tendency, Soil *soil, Grid *grid, double delta_t, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irrev, Config_info *config_info, int no_rk_step)
{
	/*
	This function manages the calculation of the explicit scalar tendencies.
	*/
	
	/*
	Firstly, some things need to prepared.
	--------------------------------------
	*/
	// declaring needed variables
    int h_index, layer_index, diff_switch;
    diff_switch = 0;
    double c_v_cond, tracer_heating, latent_heating_weight, density_total_weight;
    
    // determining the RK weights
    double old_weight[NO_OF_CONSTITUENTS];
    double new_weight[NO_OF_CONSTITUENTS];
    for (int i = 0; i < NO_OF_CONSTITUENTS; ++i)
    {
		new_weight[i] = 1;
		if (no_rk_step == 1 && i != NO_OF_CONDENSED_CONSTITUENTS)
		{
			new_weight[i] = 0.5;
		}
		old_weight[i] = 1 - new_weight[i];
    }
    
	// Temperature diffusion gets updated here, but only at the first RK step and if heat conduction is switched on.
	if (no_rk_step == 0 && (config_info -> temperature_diff_h == 1 || config_info -> temperature_diff_v == 1))
	{
	    // The diffusion of the temperature depends on its gradient.
		grad(diagnostics -> temperature_gas, diagnostics -> vector_field_placeholder, grid);
		// Now we need to calculate the temperature diffusion coefficients.
	    calc_temp_diffusion_coeffs(state, config_info, irrev, diagnostics, delta_t, grid);
		// Now the diffusive temperature flux density can be obtained.
	    scalar_times_vector_h(irrev -> scalar_diffusion_coeff_numerical_h, diagnostics -> vector_field_placeholder, diagnostics -> flux_density, grid);
	    if (config_info -> temperature_diff_v == 1)
	    {
	    	scalar_times_vector_v(irrev -> scalar_diffusion_coeff_numerical_v, diagnostics -> vector_field_placeholder, diagnostics -> flux_density, grid);
	    }
	    // The divergence of the diffusive temperature flux density is the diffusive temperature heating.
	    divv_h(diagnostics -> flux_density, irrev -> temperature_diffusion_heating, grid);
	    // the vertical divergence is only needed if the vertical temperature diffusion is switched on
	    if (config_info -> temperature_diff_v == 1)
	    {
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
        // This is the mass advection, which needs to be carried out for all constituents.
        // -------------------------------------------------------------------------------
		scalar_times_vector_h(&state -> rho[i*NO_OF_SCALARS], state -> wind, diagnostics -> flux_density, grid);
		if (i == NO_OF_CONDENSED_CONSTITUENTS)
		{
        	divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
		}
		else
		{
        	divv_h_limited(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid, &state -> rho[i*NO_OF_SCALARS], delta_t);
		}
		// horizontal mass diffusion, only for gaseous tracers
		if (config_info -> tracer_diff_h == 1 && i > NO_OF_CONDENSED_CONSTITUENTS)
		{
			diff_switch = 1;
			grad(&state -> rho[i*NO_OF_SCALARS], diagnostics -> vector_field_placeholder, grid);
			divv_h(diagnostics -> vector_field_placeholder, diagnostics -> scalar_field_placeholder, grid);
		}
		// vertical mass diffusion, only for gaseous tracers
		if (config_info -> tracer_diff_v == 1 && i > NO_OF_CONDENSED_CONSTITUENTS)
		{
			diff_switch = 1;
			if (config_info -> tracer_diff_h == 0)
			{
				grad_vert_cov(&state -> rho[i*NO_OF_SCALARS], diagnostics -> vector_field_placeholder, grid);
			}
			add_vertical_divv(diagnostics -> vector_field_placeholder, diagnostics -> scalar_field_placeholder, grid);
		}
		// adding the tendencies in all grid boxes
		#pragma omp parallel for private(layer_index, h_index)
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
			layer_index = j/NO_OF_SCALARS_H;
			h_index = j - layer_index*NO_OF_SCALARS_H;
			if (NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_scalar[h_index])
			{
				state_tendency -> rho[i*NO_OF_SCALARS + j]
				= old_weight[i]*state_tendency -> rho[i*NO_OF_SCALARS + j]
				+ new_weight[i]*(
				// the advection
				-diagnostics -> flux_density_divv[j])
				// the diffusion
				+ diff_switch*diagnostics -> scalar_field_placeholder[j];
				// the horizontal brute-force limiter
				if (state_old -> rho[i*NO_OF_SCALARS + j] + delta_t*state_tendency -> rho[i*NO_OF_SCALARS + j] < 0)
				{
					state_tendency -> rho[i*NO_OF_SCALARS + j] = -state_old -> rho[i*NO_OF_SCALARS + j]/delta_t;
				}
		    }
	    }
	    
		// explicit rho*theta integration
		// ------------------------------
		if (i == NO_OF_CONDENSED_CONSTITUENTS)
		{
			// Determining the potential temperature of the constituent at hand.
			#pragma omp parallel for
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				if (state -> rho[i*NO_OF_SCALARS + j] != 0)
				{
					diagnostics -> scalar_field_placeholder[j] = state -> rhotheta[j]/state -> rho[i*NO_OF_SCALARS + j];
				}
				else
				{
					diagnostics -> scalar_field_placeholder[j] = 0;
				}
			}
			scalar_times_vector_h(diagnostics -> scalar_field_placeholder, diagnostics -> flux_density, diagnostics -> flux_density, grid);
			divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
			// adding the tendencies in all grid boxes
			#pragma omp parallel for private(layer_index, h_index, tracer_heating, latent_heating_weight, density_total_weight)
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				layer_index = j/NO_OF_SCALARS_H;
				h_index = j - layer_index*NO_OF_SCALARS_H;
				if (NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_scalar[h_index])
				{
					// determining the heating rate that comes from the tracers
					tracer_heating = 0;
					density_total_weight = density_total(state, j);
					latent_heating_weight = 1;
					if (config_info -> assume_lte == 0)
					{
						latent_heating_weight = state -> rho[i*NO_OF_SCALARS + j]/density_gas(state, j);
						// this is not yet implemented
					}
					if (config_info -> assume_lte == 1)
					{
						latent_heating_weight = state -> rho[i*NO_OF_SCALARS + j]/density_total_weight;
						for (int k = 0; k < NO_OF_CONDENSED_CONSTITUENTS; ++k)
						{
							tracer_heating += irrev -> constituent_heat_source_rates[k*NO_OF_SCALARS + j];
						}
					}
					state_tendency -> rhotheta[j]
					= old_weight[i]*state_tendency -> rhotheta[j]
					+ new_weight[i]*(
					// the advection (resolved transport)
					-diagnostics -> flux_density_divv[j]
					// the diabatic forcings
					// weighting factor
					+ state -> rho[i*NO_OF_SCALARS + j]/density_total_weight*(
					// dissipation of molecular + turbulent momentum diffusion
					irrev -> heating_diss[j]
					// molecular + turbulent heat transport
					+ irrev -> temperature_diffusion_heating[j]
					// radiation
					+ forcings -> radiation_tendency[j]
					// this has to be divided by the c_p*exner
					)/(spec_heat_capacities_p_gas(0)*(grid -> exner_bg[j] + state -> exner_pert[j]))
					// phase transitions
					+ latent_heating_weight*tracer_heating
					/(spec_heat_capacities_p_gas(0)*(grid -> exner_bg[j] + state -> exner_pert[j])));
					// sensible heat in the lowest layer
					if (layer_index == NO_OF_LAYERS - 1 - grid -> no_of_shaded_points_scalar[h_index])
					{
						state_tendency -> rhotheta[j]
						// the minus-sign is correct (the quantity itself refers to soil)
						-= soil -> power_flux_density_sensible[h_index]/(spec_heat_capacities_p_gas(0)*(grid -> exner_bg[j] + state -> exner_pert[j]))
						/(grid -> z_vector[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER - NO_OF_SCALARS_H + h_index] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + h_index]);
					}
				 }
			}
		}
    
		// This is the integration of the "density x temperature" fields. It only needs to be done for condensed constituents.
		// -------------------------------------------------------------------------------------------------------------------
		if (i < NO_OF_CONDENSED_CONSTITUENTS && config_info -> assume_lte == 0)
		{
			// The constituent velocity has already been calculated.
		    scalar_times_vector_h(&state -> condensed_density_temperatures[i*NO_OF_SCALARS], state -> wind, diagnostics -> flux_density, grid);
		    divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
			// adding the tendencies in all grid boxes
			#pragma omp parallel for private(layer_index, h_index, c_v_cond)
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				layer_index = j/NO_OF_SCALARS_H;
				h_index = j - layer_index*NO_OF_SCALARS_H;
				if (NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_scalar[h_index])
				{
					c_v_cond = ret_c_v_cond(i, 0, state -> condensed_density_temperatures[i*NO_OF_SCALARS + j]/(EPSILON_SECURITY + state -> rho[i*NO_OF_SCALARS + j]));
					state_tendency -> condensed_density_temperatures[i*NO_OF_SCALARS + j]
					= old_weight[i]*state_tendency -> condensed_density_temperatures[i*NO_OF_SCALARS + j]
					+ new_weight[i]*(
					// the advection
					-diagnostics -> flux_density_divv[j]
					// the source terms
					+ state -> rho[i*NO_OF_SCALARS + j]/(EPSILON_SECURITY + c_v_cond*density_total(state, j))
					*(irrev -> temperature_diffusion_heating[j] + irrev -> heating_diss[j] + forcings -> radiation_tendency[j])
					+ 1/c_v_cond*irrev -> constituent_heat_source_rates[i*NO_OF_SCALARS + j]
					+ state -> condensed_density_temperatures[i*NO_OF_SCALARS + j]*(irrev -> mass_source_rates[i*NO_OF_SCALARS + j]));
				}
			}
		}
	} // constituent loop
	
	return 0;
}

int moisturizer(State *state, double delta_t, Diagnostics *diagnostics, Irreversible_quantities *irrev, Config_info *config_info, Grid *grid)
{
	/*
	This function manages the calculation of the phase transition rates.
	*/
	
	// Only if we have multiple constituents, moisture needs to be included.
	if (NO_OF_CONSTITUENTS > 1)
	{
		// calculating the source rates
	    calc_h2otracers_source_rates(
	    irrev -> mass_source_rates,
	    irrev -> constituent_heat_source_rates,
	    state -> rho,
	    state -> condensed_density_temperatures,
	    diagnostics -> temperature_gas,
	    NO_OF_SCALARS,
	    2*delta_t,
	    config_info -> assume_lte,
	    grid -> is_land,
	    NO_OF_LAYERS);
	    int layer_index, h_index;
	    // loop over all constituents
		for (int i = 0; i < NO_OF_CONSTITUENTS; ++i)
		{
			// the main gaseous constituent has no source rates
			if (i != NO_OF_CONDENSED_CONSTITUENTS)
			{
				#pragma omp parallel for private(layer_index, h_index)
				for (int j = 0; j < NO_OF_SCALARS; ++j)
				{
					layer_index = j/NO_OF_SCALARS_H;
					h_index = j - layer_index*NO_OF_SCALARS_H;
					// check for shading
					if (NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_scalar[h_index])
					{
						if (i < NO_OF_CONDENSED_CONSTITUENTS)
						{
							state -> rho[i*NO_OF_SCALARS + j] = state -> rho[i*NO_OF_SCALARS + j] + delta_t*irrev -> mass_source_rates[i*NO_OF_SCALARS + j];
						}
						// for the gaseous constituents (apart from the main one), an index shift is necessary
						else
						{
							state -> rho[i*NO_OF_SCALARS + j] = state -> rho[i*NO_OF_SCALARS + j] + delta_t*irrev -> mass_source_rates[(i - 1)*NO_OF_SCALARS + j];
						}
					}
				}
			}
		}
	}
	return 0;	
}






