/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains functions calculating everything related to phase transition rates.
*/

#include <math.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../subgrid_scale/subgrid_scale.h"
#include "constituents.h"

int calc_h2otracers_source_rates(State *state, Diagnostics *diagnostics, Grid *grid, Config *config, Irreversible_quantities *irrev, double delta_t)
{
	/*
	This function calculates phase transition rates and associated heat source rates.
	It assumes the following order for the constituents:
	precipitating ice - precipitating liquid water - cloud ice - liquid cloud water - dry air - water vapour
	*/
	
    double diff_density, phase_trans_density, saturation_pressure, water_vapour_pressure, solid_temperature, liquid_temperature,
    layer_thickness, diff_density_sfc, saturation_pressure_sfc, dry_pressure, air_pressure;
    
    //  maximum cloud water content in (kg cloud)/(kg dry air).
    double maximum_cloud_water_content = 0.2e-3;
    
    // loop over all grid boxes
    int layer_index, h_index;
    #pragma omp parallel for private(diff_density, phase_trans_density, saturation_pressure, water_vapour_pressure, solid_temperature, liquid_temperature, layer_index, h_index, layer_thickness, diff_density_sfc, saturation_pressure_sfc, dry_pressure, air_pressure)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	/*
    	Preparation
    	-----------
    	*/
    	layer_index = i/NO_OF_SCALARS_H;
    	// determining the temperature of the cloud ice
    	if (state -> rho[2*NO_OF_SCALARS + i] < EPSILON_SECURITY)
    	{
    		solid_temperature = T_0;
		}
		else
		{
			solid_temperature = diagnostics -> temperature_gas[i];
		}
		
		// determining the temperature of the liquid cloud water
    	if (state -> rho[3*NO_OF_SCALARS + i] < EPSILON_SECURITY)
    	{
    		liquid_temperature = T_0;
		}
		else
		{
			liquid_temperature = diagnostics -> temperature_gas[i];
		}
		
		// determining the saturation pressure
		// "positive" temperatures (the saturation pressure is different over water compared to over ice)
        if (diagnostics -> temperature_gas[i] >= T_0)
    	{
            saturation_pressure = saturation_pressure_over_water(diagnostics -> temperature_gas[i]);
		}
		// "negative" temperatures
        else
    	{
            saturation_pressure = saturation_pressure_over_ice(diagnostics -> temperature_gas[i]);
		}
		
		// determining the water vapour pressure (using the EOS)
        water_vapour_pressure = state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i]*specific_gas_constants(1)*diagnostics -> temperature_gas[i];
		
		// determining the water vapour pressure (using the EOS)
        dry_pressure = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*specific_gas_constants(0)*diagnostics -> temperature_gas[i];
        
        // calculating the total air pressure
        air_pressure = dry_pressure + water_vapour_pressure;
        
        // multiplying the saturation pressure by the enhancement factor
        if (diagnostics -> temperature_gas[i] >= T_0)
    	{
            saturation_pressure = enhancement_factor_over_water(air_pressure)*saturation_pressure;
		}
		// "negative" temperatures
        else
    	{
            saturation_pressure = enhancement_factor_over_ice(air_pressure)*saturation_pressure;
		}
        
    	// the amount of water vapour that the air can still take 
        diff_density = (saturation_pressure - water_vapour_pressure)/(specific_gas_constants(1)*diagnostics -> temperature_gas[i]);
        
    	/*
    	Clouds
    	------
    	*/
        // the case where the air is not over-saturated
        if (diff_density >= 0.0)
        {
            // temperature >= 0 °C
            if (diagnostics -> temperature_gas[i] >= T_0)
            {
            	// It is assumed that the still present ice vanishes within one time step.
                irrev -> mass_source_rates[2*NO_OF_SCALARS + i] = -state -> rho[2*NO_OF_SCALARS + i]/delta_t;
                
                /*
                The amount of liquid water per volume that will evaporate.
                In case the air cannot take all the water, not everything will evaporate.
                */
                phase_trans_density = fmin(state -> rho[3*NO_OF_SCALARS + i], diff_density);
                
                /*
                The source rate for the liquid water consists of two terms:
                1.) the evaporation
                 2.) the melting of ice
                */
                irrev -> mass_source_rates[3*NO_OF_SCALARS + i] = (state -> rho[2*NO_OF_SCALARS + i] - phase_trans_density)/delta_t;
                
                // the tendency for the water vapour
                irrev -> mass_source_rates[4*NO_OF_SCALARS + i] = phase_trans_density/delta_t;
                
                // the heat source rates acting on the ice
                irrev -> phase_trans_heating_rate[i] = irrev -> mass_source_rates[2*NO_OF_SCALARS + i]*phase_trans_heat(2, solid_temperature);
                
                // the heat source rates acting on the liquid water
                irrev -> phase_trans_heating_rate[i] +=
                // the evaporation
                -phase_trans_density*phase_trans_heat(0, T_0)/delta_t;
            }
            // temperature < 0 °C
            else
            {
            	// Everything that can sublimate will sublimate.
                phase_trans_density = fmin(state -> rho[2*NO_OF_SCALARS + i], diff_density);
                
                /*
                the tendency for the ice contains two terms:
                1.) the freezing
                2.) the phase transition through sublimation
                */
                irrev -> mass_source_rates[2*NO_OF_SCALARS + i] = (state -> rho[3*NO_OF_SCALARS + i] - phase_trans_density)/delta_t;
                
            	// It is assumed that the still present liquid water vanishes within one time step.
                irrev -> mass_source_rates[3*NO_OF_SCALARS + i] = -state -> rho[3*NO_OF_SCALARS + i]/delta_t;
                
                // the tendency for the water vapour
                irrev -> mass_source_rates[4*NO_OF_SCALARS + i] = phase_trans_density/delta_t;
                
                // the heat source rates acting on the ice
                irrev -> phase_trans_heating_rate[i] = (
                // the freezing
                state -> rho[3*NO_OF_SCALARS + i]*phase_trans_heat(2, solid_temperature)
                // the sublimation
                - phase_trans_density*phase_trans_heat(1, solid_temperature))/delta_t;
            }
        }
        // the case where the air is over-saturated
        else
        {
        	// the vanishing of water vapour through the phase transition
            irrev -> mass_source_rates[4*NO_OF_SCALARS + i] = diff_density/delta_t;
            // temperature >= 0 °C
            if (diagnostics -> temperature_gas[i] >= T_0)
            {
            	// It is assumed that the still present ice vanishes within one time step.
                irrev -> mass_source_rates[2*NO_OF_SCALARS + i] = -state -> rho[2*NO_OF_SCALARS + i]/delta_t;
                
                /*
                The source rate for the liquid water consists of two terms:
                1.) the condensation
                2.) the melting of ice
				*/
				irrev -> mass_source_rates[3*NO_OF_SCALARS + i] = (-diff_density + state -> rho[2*NO_OF_SCALARS + i])/delta_t;
                
                // the heat source rates acting on the ice
                irrev -> phase_trans_heating_rate[i] = irrev -> mass_source_rates[2*NO_OF_SCALARS + i]*phase_trans_heat(2, solid_temperature);
                
                // the heat source rates acting on the liquid water
                irrev -> phase_trans_heating_rate[i] +=
                // it is only affected by the condensation
                -diff_density*phase_trans_heat(0, liquid_temperature)/delta_t;
            }
            // temperature < 0 °C
            else
            {
            	/*
                The source rate for the cloud ice consists of two terms:
                1.) the resublimation
                2.) the freezing of cloud water
                */
                irrev -> mass_source_rates[2*NO_OF_SCALARS + i] = (-diff_density + state -> rho[3*NO_OF_SCALARS + i])/delta_t;
                
                // It is assumed that the liquid water disappears within one time step.
                irrev -> mass_source_rates[3*NO_OF_SCALARS + i] = -state -> rho[3*NO_OF_SCALARS + i]/delta_t;
                
                // the heat source rates acting on the ice
                irrev -> phase_trans_heating_rate[i] =
                // the component through the resublimation
                (-diff_density*phase_trans_heat(1, solid_temperature)
                // the component through freezing
                + state -> rho[3*NO_OF_SCALARS + i]*phase_trans_heat(2, solid_temperature))/delta_t;
            }
        }
        
    	/*
    	Precipitation
    	-------------
    	*/
        irrev -> mass_source_rates[i] = 0.0;
        irrev -> mass_source_rates[NO_OF_SCALARS + i] = 0.0;
        // snow
        // this only happens if the air is saturated
        if (diagnostics -> temperature_gas[i] < T_0 && diff_density <= 0.0)
        {
        	irrev -> mass_source_rates[i] = fmax(state -> rho[2*NO_OF_SCALARS + i] - maximum_cloud_water_content*state -> rho[4*NO_OF_SCALARS + i], 0.0)/delta_t;
        	// the snow creation comes at the cost of cloud ice particles
        	irrev -> mass_source_rates[2*NO_OF_SCALARS + i] -= irrev -> mass_source_rates[i];
        }
    	// rain
        // this only happens if the air is saturated
        else if (diagnostics -> temperature_gas[i] >= T_0 && diff_density <= 0.0)
        {
        	irrev -> mass_source_rates[NO_OF_SCALARS + i] = fmax(state -> rho[3*NO_OF_SCALARS + i] - maximum_cloud_water_content*state -> rho[4*NO_OF_SCALARS + i], 0.0)/delta_t;
        	// the rain creation comes at the cost of cloud water particles
        	irrev -> mass_source_rates[3*NO_OF_SCALARS + i] -= irrev -> mass_source_rates[NO_OF_SCALARS + i];
        }
        
        // turning of snow to rain
        if (diagnostics -> temperature_gas[i] >= T_0 && state -> rho[i] > 0.0)
        {
        	irrev -> mass_source_rates[i] = -state -> rho[i]/delta_t;
        	irrev -> mass_source_rates[NO_OF_SCALARS + i] -= irrev -> mass_source_rates[i];
        	irrev -> phase_trans_heating_rate[i] += irrev -> mass_source_rates[i]*phase_trans_heat(2, T_0);
        }
        // turning of rain to snow
        if (diagnostics -> temperature_gas[i] < T_0 && state -> rho[NO_OF_SCALARS + i] > 0.0)
        {
        	irrev -> mass_source_rates[NO_OF_SCALARS + i] = -state -> rho[NO_OF_SCALARS + i]/delta_t;
        	irrev -> mass_source_rates[i] -= irrev -> mass_source_rates[NO_OF_SCALARS + i];
        	irrev -> phase_trans_heating_rate[i] += -irrev -> mass_source_rates[NO_OF_SCALARS + i]*phase_trans_heat(2, T_0);
        }
        
        /*
        Surface effects
        ---------------
        */
        if (layer_index == NO_OF_LAYERS - 1 && config -> sfc_phase_trans == 1)
        {
	    	h_index = i - layer_index*NO_OF_SCALARS_H;
	    	
	    	// evaporation and latent heat rates
    		if (grid -> is_land[h_index] == 0)
        	{
        		// saturation pressure at surface temperature
        		if (state -> temperature_soil[h_index] >= T_0)
        		{
        			saturation_pressure_sfc = saturation_pressure_over_water(state -> temperature_soil[h_index]);
		    		saturation_pressure_sfc = enhancement_factor_over_water(air_pressure)*saturation_pressure_sfc;
        		}
        		else
        		{
        			saturation_pressure_sfc = saturation_pressure_over_ice(state -> temperature_soil[h_index]);
		    		saturation_pressure_sfc = enhancement_factor_over_ice(air_pressure)*saturation_pressure_sfc;
        		}
        		
        		// difference water vapour density between saturation at ground temperature and actual absolute humidity in the lowest model layer
        		diff_density_sfc = saturation_pressure_sfc/(specific_gas_constants(1)*state -> temperature_soil[h_index])
        		- state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i];
        		
        		// the thickness of the lowest model layer (we need it as a result of Guass' theorem)
        		layer_thickness = grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + h_index] - grid -> z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + h_index];
        		
        		// evporation, sublimation
		    	irrev -> mass_source_rates[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i] += fmax(0.0, diff_density_sfc/diagnostics -> scalar_flux_resistance[h_index])/layer_thickness;
		    	
		    	// calculating the latent heat flux density affecting the surface
        		if (state -> temperature_soil[h_index] >= T_0)
        		{
        			diagnostics -> power_flux_density_latent[h_index] = -phase_trans_heat(0, state -> temperature_soil[h_index])
        			*fmax(0.0, diff_density_sfc/diagnostics -> scalar_flux_resistance[h_index]);
        		}
        		else
        		{
        			diagnostics -> power_flux_density_latent[h_index] = -phase_trans_heat(1, state -> temperature_soil[h_index])
        			*fmax(0.0, diff_density_sfc/diagnostics -> scalar_flux_resistance[h_index]);
        		}
        	}
        }
    }
	
    return 0;
}














