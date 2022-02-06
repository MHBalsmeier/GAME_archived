/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains functions calculating everything related to phase transition rates.
*/

#include <stdlib.h>
#include <math.h>
#include "tracers.h"

int calc_h2otracers_source_rates(double mass_source_rates[], double heat_source_rates[], double densities[],
double tracer_temperature_densities[], double temperature[], int number_of_scalars, double delta_t, int assume_lte, int is_land[], int number_of_layers,
double z_scalar[], double z_vector[], double v_squared[], int no_of_vector_per_layer, double roughness_length[], double temperature_soil[],
double power_flux_density_sensible[], double power_flux_density_latent[], int soil_on)
{
	/*
	This function calculates phase transition rates and associated heat source rates.
	It assumes the following order for the constituents:
	precipitating ice - precipitating liquid water - cloud ice - liquid cloud water - dry air - water vapour
	*/
    double diff_density, phase_trans_density, saturation_pressure, water_vapour_pressure, solid_temperature, liquid_temperature,
    flux_resistance, layer_thickness, diff_density_sfc, saturation_pressure_sfc;
    
    //  maximum cloud water content in (kg cloud)/(kg dry air).
    double maximum_cloud_water_content = 0.2e-3;
    
    int number_of_scalars_h = number_of_scalars/number_of_layers;
    
    // loop over all grid boxes
    int layer_index, h_index;
    #pragma omp parallel for private(diff_density, phase_trans_density, saturation_pressure, water_vapour_pressure, solid_temperature, liquid_temperature, layer_index, h_index, flux_resistance, layer_thickness, diff_density_sfc, saturation_pressure_sfc)
    for (int i = 0; i < number_of_scalars; ++i)
    {
    	layer_index = i/number_of_scalars_h;
    	// determining the temperature of the cloud ice
    	if (densities[2*number_of_scalars + i] < EPSILON_SECURITY)
    	{
    		solid_temperature = T_0;
		}
		else if (assume_lte == 1)
		{
			solid_temperature = temperature[i];
		}
    	else
    	{
    		solid_temperature = tracer_temperature_densities[2*number_of_scalars + i]/densities[2*number_of_scalars + i];
		}
		
		// determining the temperature of the liquid cloud water
    	if (densities[3*number_of_scalars + i] < EPSILON_SECURITY)
    	{
    		liquid_temperature = T_0;
		}
		else if (assume_lte == 1)
		{
			liquid_temperature = temperature[i];
		}
    	else
    	{
    		liquid_temperature = tracer_temperature_densities[3*number_of_scalars + i]/densities[3*number_of_scalars + i];
		}
		
		// determining the saturation pressure
		// "positive" temperatures (the saturation pressure is different over water compared to over ice)
        if (temperature[i] >= T_0)
    	{
            saturation_pressure = saturation_pressure_over_water(temperature[i]);
		}
		// "negative" temperatures
        else
    	{
            saturation_pressure = saturation_pressure_over_ice(temperature[i]);
		}
		
		// determining the water vapour pressure (using the EOS)
        water_vapour_pressure = densities[5*number_of_scalars + i]*specific_gas_constants_lookup(1)*temperature[i];
        
    	// the amount of water vapour that the air can still take 
        diff_density = (saturation_pressure - water_vapour_pressure)/(specific_gas_constants_lookup(1)*temperature[i]);
        
        // the case where the air is not over-saturated
        if (diff_density >= 0)
        {
            // temperature >= 0 °C
            if (temperature[i] >= T_0)
            {
            	// It is assumed that the still present ice vanishes within one time step.
                mass_source_rates[2*number_of_scalars + i] = -densities[2*number_of_scalars + i]/delta_t;
                
                /*
                The amount of liquid water per volume that will evaporate.
                In case the air cannot take all the water, not everything will evaporate.
                */
                phase_trans_density = fmin(densities[3*number_of_scalars + i], diff_density);
                
                /*
                The source rate for the liquid water consists of two terms:
                1.) the evaporation
                 2.) the melting of ice
                */
                mass_source_rates[3*number_of_scalars + i] = (densities[2*number_of_scalars + i] - phase_trans_density)/delta_t;
                
                // the tendency for the water vapour
                mass_source_rates[4*number_of_scalars + i] = phase_trans_density/delta_t;
                
                // the heat source rates acting on the ice
                heat_source_rates[2*number_of_scalars + i] = mass_source_rates[2*number_of_scalars + i]*phase_trans_heat(2, solid_temperature);
                
                // the heat source rates acting on the liquid water
                heat_source_rates[3*number_of_scalars + i] =
                // the evaporation
                -phase_trans_density*phase_trans_heat(0, T_0)/delta_t;
            }
            // temperature < 0 °C
            else
            {
            	// Everything that can sublimate will sublimate.
                phase_trans_density = fmin(densities[2*number_of_scalars + i], diff_density);
                
                /*
                the tendency for the ice contains two terms:
                1.) the freezing
                2.) the phase transition through sublimation
                */
                mass_source_rates[2*number_of_scalars + i] = (densities[3*number_of_scalars + i] - phase_trans_density)/delta_t;
                
            	// It is assumed that the still present liquid water vanishes within one time step.
                mass_source_rates[3*number_of_scalars + i] = -densities[3*number_of_scalars + i]/delta_t;
                
                // the tendency for the water vapour
                mass_source_rates[4*number_of_scalars + i] = phase_trans_density/delta_t;
                
                // the heat source rates acting on the ice
                heat_source_rates[2*number_of_scalars + i] = (
                // the freezing
                densities[3*number_of_scalars + i]*phase_trans_heat(2, solid_temperature)
                // the sublimation
                - phase_trans_density*phase_trans_heat(1, solid_temperature))/delta_t;
                
                // the heat source rates acting on the liquid water
                heat_source_rates[3*number_of_scalars + i] = 0;
            }
        }
        // the case where the air is over-saturated
        else
        {
        	// the vanishing of water vapour through the phase transition
            mass_source_rates[4*number_of_scalars + i] = diff_density/delta_t;
            // temperature >= 0 °C
            if (temperature[i] >= T_0)
            {
            	// It is assumed that the still present ice vanishes within one time step.
                mass_source_rates[2*number_of_scalars + i] = -densities[2*number_of_scalars + i]/delta_t;
                
                /*
                The source rate for the liquid water consists of two terms:
                1.) the condensation
                2.) the melting of ice
				*/
				mass_source_rates[3*number_of_scalars + i] = (-diff_density + densities[2*number_of_scalars + i])/delta_t;
                
                // the heat source rates acting on the ice
                heat_source_rates[2*number_of_scalars + i] = -densities[2*number_of_scalars + i]*phase_trans_heat(2, solid_temperature)/delta_t;
                
                // the heat source rates acting on the liquid water
                heat_source_rates[3*number_of_scalars + i] =
                // it is only affected by the condensation
                -diff_density*phase_trans_heat(0, liquid_temperature)/delta_t;
            }
            // temperature < 0 °C
            else
            {
            	/*
                The source rate for the ice consists of two terms:
                1.) the resublimation
                2.) the melting of ice
                */
                mass_source_rates[2*number_of_scalars + i] = (-diff_density + densities[3*number_of_scalars + i])/delta_t;
                
                // It is assumed that the liquid water disappears within one time step.
                mass_source_rates[3*number_of_scalars + i] = -densities[3*number_of_scalars + i]/delta_t;
                
                // the heat source rates acting on the ice
                heat_source_rates[2*number_of_scalars + i] =
                // the component through the resublimation
                (-diff_density*phase_trans_heat(1, solid_temperature)
                // the component through freezing
                + densities[3*number_of_scalars + i]*phase_trans_heat(2, solid_temperature))/delta_t;
                
                // the heat source rates acting on the liquid water
                heat_source_rates[3*number_of_scalars + i] = 0;
            }
        }
        
        // creation of precipitation
        // snow
        mass_source_rates[i] = fmax(densities[2*number_of_scalars + i] - maximum_cloud_water_content*densities[4*number_of_scalars + i], 0)/delta_t;
        // the snow creation comes at the cost of cloud ice particles
        mass_source_rates[2*number_of_scalars + i] -= mass_source_rates[i];
        // rain
        mass_source_rates[number_of_scalars + i] = fmax(densities[3*number_of_scalars + i] - maximum_cloud_water_content*densities[4*number_of_scalars + i], 0)/delta_t;
        // the rain creation comes at the cost of cloud water particles
        mass_source_rates[3*number_of_scalars + i] -= mass_source_rates[number_of_scalars + i];
        
        // turning of snow to rain
        if (temperature[i] > T_0 && densities[i] > 0)
        {
        	mass_source_rates[i] = -densities[i]/delta_t;
        	mass_source_rates[number_of_scalars + i] -= mass_source_rates[i];
        }
        
        // surface effects
        if (layer_index == number_of_layers - 1 && soil_on == 1)
        {
	    	h_index = i - layer_index*number_of_scalars_h;	    	
	    	
    		// flux resistance
	    	flux_resistance = sfc_flux_resistance(pow(v_squared[i], 0.5), z_scalar[i] - z_vector[(layer_index + 1)*no_of_vector_per_layer + h_index], roughness_length[h_index]);
			// sensible heat flux density through the surface (towards the surface)
			power_flux_density_sensible[h_index] = densities[4*number_of_scalars + i]
			*spec_heat_capacities_v_gas_lookup(0)*(temperature[i] - temperature_soil[h_index])/flux_resistance;
	    	
	    	// evaporation and latent heat rates
    		if (is_land[h_index] == 0)
        	{
        		// saturation pressure at surface temperature
        		if (temperature_soil[h_index] >= T_0)
        		{
        			saturation_pressure_sfc = saturation_pressure_over_water(temperature_soil[h_index]);
        		}
        		else
        		{
        			saturation_pressure_sfc = saturation_pressure_over_ice(temperature_soil[h_index]);
        		}
        		// difference water vapour density between saturation at ground temperature and actual absolute humidity in the lowest model layer
        		diff_density_sfc = saturation_pressure_sfc/(specific_gas_constants_lookup(1)*temperature_soil[h_index])
        		- water_vapour_pressure/(specific_gas_constants_lookup(1)*temperature[i]);
        		// the thickness of the lowest model layer (we need it as a result of Guass' theorem)
        		layer_thickness = z_vector[layer_index*no_of_vector_per_layer + h_index] - z_vector[(layer_index + 1)*no_of_vector_per_layer + h_index];
		    	mass_source_rates[4*number_of_scalars + i] += fmax(0, diff_density_sfc/flux_resistance)/layer_thickness;
		    	// calculating the latent heat flux density affecting the surface
        		if (temperature_soil[h_index] >= T_0)
        		{
        			power_flux_density_latent[h_index] = -phase_trans_heat(0, temperature_soil[h_index])*fmax(0, diff_density_sfc/flux_resistance);
        		}
        		else
        		{
        			power_flux_density_latent[h_index] = -phase_trans_heat(1, temperature_soil[h_index])*fmax(0, diff_density_sfc/flux_resistance);
        		}
        	}
        }
    }
    return 0;
}

double phase_trans_heat(int direction, double temperature)
{
    /*
    directions:
    0:  gas to liquid
    1:  gas to solid
    2: liquid to solid
    */
    double result;
    if (direction == 0)
    {
        result = 2257000;
	}
    if (direction == 1)
    {
        result = 2257000 + 333500;
	}
    if (direction == 2)
    {
        result = 333500;
	}
    return result;
}

double sink_velocity(int solid_or_liquid, double radius, double air_density)
{
	double dry_air_kinematic_viscosity = 14.8e-6;
	double reynolds_crit = 10;
	double drag_coeff = 1;
	
	// First of all, a laminar sink velocity is calculated from the Stokes law.
	double laminar_velocity_candidate = 0;
	// The solid case.
	if (solid_or_liquid == 0)
	{
		laminar_velocity_candidate = 2*M_PI*pow(radius, 2)*DENSITY_WATER*GRAVITY/(9*M_PI*air_density*dry_air_kinematic_viscosity);
	}
	
	// The liquid case.
	if (solid_or_liquid == 1)
	{
		laminar_velocity_candidate = 2*M_PI*pow(radius, 2)*DENSITY_WATER*GRAVITY/(9*M_PI*air_density*dry_air_kinematic_viscosity);
	}
	
	// calculating the Reynolds number resulting from the laminar velocity
	double reynolds_from_laminar;
	reynolds_from_laminar = laminar_velocity_candidate*radius/dry_air_kinematic_viscosity;
	
	// calculating the resulting sink velocity
	double result;
	// the laminar case
	if (reynolds_from_laminar <= reynolds_crit)
	{
		result = laminar_velocity_candidate;
	}
	// the turbulent case
	else
	{
		result = pow(8*radius*DENSITY_WATER*GRAVITY/(3*air_density*drag_coeff), 0.5);
	}
	
    return result;
}


