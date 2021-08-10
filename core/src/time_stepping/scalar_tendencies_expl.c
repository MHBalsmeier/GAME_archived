/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/game
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

int scalar_tendencies_expl(State *state, State *state_tendency, Soil *soil, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irrev, Config_info *config_info, int no_rk_step, Vector_field wind_advect_tracer)
{
	/*
	Firstly, some things need to prepared.
	--------------------------------------
	*/
	// declaring needed variables
    int h_index, layer_index;
    double c_v_cond;
    
    // determining the weights for the RK stepping
    double old_weight, new_weight, tracer_heating, density_gas_weight, density_total_weight;
    new_weight = 1;
    if (no_rk_step == 1)
    {
    	new_weight = 1;
    }
	old_weight = 1 - new_weight;
    
	// phase transitions are only updated at the first RK step
	if (NO_OF_CONSTITUENTS == 4)
	{
	    calc_h2otracers_source_rates(
	    irrev -> constituent_mass_source_rates,
	    irrev -> constituent_heat_source_rates,
	    state -> rho,
	    state -> condensed_density_temperatures,
	    diagnostics -> temperature_gas,
	    NO_OF_SCALARS,
	    delta_t,
	    config_info -> assume_lte);
	}
	
	// Temperature diffusion gets updated here, but only at the first RK step and if heat conduction is switched on.
	if (no_rk_step == 0 && (config_info -> temperature_diff_h == 1 || config_info -> temperature_diff_v == 1))
	{
	    // The diffusion of the temperature depends on its gradient.
		grad(diagnostics -> temperature_gas, diagnostics -> vector_field_placeholder, grid);
		// Now we need to calculate the temperature diffusion coefficients.
	    calc_temp_diffusion_coeffs(state, config_info, irrev, diagnostics, delta_t, grid);
		// Now the diffusive temperature flux density can be obtained.
	    scalar_times_vector_scalar_h(irrev -> scalar_diffusion_coeff_numerical_h, diagnostics -> vector_field_placeholder, diagnostics -> flux_density, grid);
	    if (config_info -> temperature_diff_v == 1)
	    {
	    	scalar_times_vector_scalar_v(irrev -> scalar_diffusion_coeff_numerical_v, diagnostics -> vector_field_placeholder, diagnostics -> flux_density, grid);
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
		// Separating the mass density of the constituent at hand.
		#pragma omp parallel for
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
		    diagnostics -> scalar_field_placeholder[j] = state -> rho[i*NO_OF_SCALARS + j];
	    }
        
        // This is the mass advection, which needs to be carried out for all constituents.
        // -------------------------------------------------------------------------------
		scalar_times_vector(diagnostics -> scalar_field_placeholder, state -> wind, diagnostics -> flux_density, grid);
        divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
		// adding the tendencies in all grid boxes
		#pragma omp parallel for private(layer_index, h_index)
		for (int j = 0; j < NO_OF_SCALARS; ++j)
		{
			layer_index = j/NO_OF_SCALARS_H;
			h_index = j - layer_index*NO_OF_SCALARS_H;
			if (NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_scalar[h_index])
			{
				state_tendency -> rho[i*NO_OF_SCALARS + j]
				= old_weight*state_tendency -> rho[i*NO_OF_SCALARS + j]
				+ new_weight*(
				// the advection
				-diagnostics -> flux_density_divv[j]
				// the phase transition rates
				+ irrev -> constituent_mass_source_rates[i*NO_OF_SCALARS + j]);
		    }
	    }
	    
		// explicit entropy integrations
		// -----------------------------
		if (i == NO_OF_CONDENSED_CONSTITUENTS)
		{
			// Determining the specific entropy of the constituent at hand.
			#pragma omp parallel for
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				if (state -> rho[i*NO_OF_SCALARS + j] != 0)
				{
					diagnostics -> scalar_field_placeholder[j] = grid -> theta_bg[j] + state -> theta_pert[j];
				}
				else
				{
					diagnostics -> scalar_field_placeholder[j] = 0;
				}
			}
			scalar_times_vector(diagnostics -> scalar_field_placeholder, diagnostics -> flux_density, diagnostics -> flux_density, grid);
			divv_h(diagnostics -> flux_density, diagnostics -> flux_density_divv, grid);
			// adding the tendencies in all grid boxes
			#pragma omp parallel for private(layer_index, h_index, tracer_heating, density_gas_weight, density_total_weight)
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				layer_index = j/NO_OF_SCALARS_H;
				h_index = j - layer_index*NO_OF_SCALARS_H;
				if (NO_OF_LAYERS - 1 - layer_index >= grid -> no_of_shaded_points_scalar[h_index])
				{
					// determining the heating rate that comes from the tracers
					tracer_heating = 0;
					density_gas_weight = 0;
					density_total_weight = 0;
					if (config_info -> assume_lte == 0)
					{
						density_gas_weight = density_gas(state, j);
						density_total_weight = density_total(state, j);
						tracer_heating = irrev -> constituent_heat_source_rates[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j];
					}
					if (config_info -> assume_lte == 1)
					{
						density_gas_weight = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j];
						density_total_weight = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + j];
						for (int k = 0; k < NO_OF_CONDENSED_CONSTITUENTS + 1; ++k)
						{
							tracer_heating += irrev -> constituent_heat_source_rates[k*NO_OF_SCALARS + j];
						}
					}
					state_tendency -> rhotheta[j]
					= old_weight*state_tendency -> rhotheta[j]
					+ new_weight*(
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
					+ radiation_tendency[j]
					// this has to be divided by the temperature (we ware in the entropy equation)
					)/diagnostics -> temperature_gas[j]
					// phase transitions
					+ tracer_heating*state -> rho[i*NO_OF_SCALARS + j]/density_gas_weight
					/diagnostics -> temperature_gas[j]);
					// sensible heat in the lowest layer
					if (layer_index == NO_OF_LAYERS - 1 - grid -> no_of_shaded_points_scalar[h_index])
					{
						state_tendency -> rhotheta[j]
						// the minus-sign is correct (the quantity itself refers to soil)
						-= new_weight*soil -> power_flux_density_sensible[j + NO_OF_SCALARS_H - NO_OF_SCALARS]/diagnostics -> temperature_gas[j]
						/(grid -> z_vector[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER - NO_OF_SCALARS_H + h_index] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + h_index]);
					}
				 }
			}
		}
    
		// This is the integration of the "density x temperature" fields. It only needs to be done for condensed constituents.
		// -------------------------------------------------------------------------------------------------------------------
		if (i < NO_OF_CONDENSED_CONSTITUENTS && config_info -> assume_lte == 0)
		{
			#pragma omp parallel for
			for (int j = 0; j < NO_OF_SCALARS; ++j)
			{
				diagnostics -> scalar_field_placeholder[j] = state -> condensed_density_temperatures[i*NO_OF_SCALARS + j];
			}
			// The constituent velocity has already been calculated.
		    scalar_times_vector(diagnostics -> scalar_field_placeholder, state -> wind, diagnostics -> flux_density, grid);
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
					= old_weight*state_tendency -> condensed_density_temperatures[i*NO_OF_SCALARS + j]
					+ new_weight*(
					// the advection
					-diagnostics -> flux_density_divv[j]
					// the source terms
					+ state -> rho[i*NO_OF_SCALARS + j]/(EPSILON_SECURITY + c_v_cond*density_total(state, j))
					*(irrev -> temperature_diffusion_heating[j] + irrev -> heating_diss[j] + radiation_tendency[j])
					+ 1/c_v_cond*irrev -> constituent_heat_source_rates[i*NO_OF_SCALARS + j]
					+ diagnostics -> scalar_field_placeholder[j]*(irrev -> constituent_mass_source_rates[i*NO_OF_SCALARS + j]));
				}
			}
		}
	} // constituent loop
	
	return 0;
}









