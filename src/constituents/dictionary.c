/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains look-up functions for properties of the atmosphere.
*/

#include <math.h>
#include "../game_constants.h"
#include "../game_types.h"
#include "constituents.h"

/*
Gas quantities
--------------
*/

/*
gaseous constituents IDs:
0: dry air
1: H2O
2: N2
3: O2
4: Ar
5: CO2
6: Ne
7: He
8: CH4
9: CO
10: O3
11: N2O
*/

double enthalpy_evaporation(double);
double enthalpy_melting(double);

double mean_particle_masses_gas(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 0.004810e-23;
	}
	if (gas_constituent_id == 1)
	{
		result = 0.002991e-23;
	}
	return result;
}

double spec_heat_capacities_v_gas(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 717.942189;
	}
	if (gas_constituent_id == 1)
	{
		result = 1396.475121;
	}
	return result;
}

double spec_heat_capacities_p_gas(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 1005.0;
	}
	if (gas_constituent_id == 1)
	{
		result = 1858.0;
	}
	return result;
}

double specific_gas_constants(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 287.057811;
	}
	if (gas_constituent_id == 1)
	{
		result = 461.524879;
	}
	return result;
}

// This follows Zdunkowski & Bott: Thermodynamics of the Atmosphere (2004), pp. 120ff.
double molar_fraction_in_dry_air(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 2)
	{
		result = 0.7809;
	}
	if (gas_constituent_id == 3)
	{
		result = 0.2095;
	}
	if (gas_constituent_id == 4)
	{
		result = 0.0093;
	}
	if (gas_constituent_id == 5)
	{
		result = 0.0003;
	}
	if (gas_constituent_id == 6)
	{
		result = 1.8e-5;
	}
	if (gas_constituent_id == 7)
	{
		result = 5.2e-6;
	}
	if (gas_constituent_id == 8)
	{
		result = 1.5e-6;
	}
	if (gas_constituent_id == 9)
	{
		result = 1.0e-7;
	}
	if (gas_constituent_id == 10)
	{
		result = 1e-6;
	}
	if (gas_constituent_id == 11)
	{
	    // https://www.epa.gov/climate-indicators/climate-change-indicators-atmospheric-concentrations-greenhouse-gases
		result = 0.3e-6;
	}
	return result;
}

double calc_o3_vmr(double z_height)
{
	/*
	calculates the ozone VMR as a function of height,
	assumes a Gaussian distribution
	*/
	
	double fwhm = 20e3;
	double z_max = 34e3;
	double max_vmr = 8.5e-6;
	
	// calculation of the result
	double sigma = fwhm/pow(8.0*log(2.0), 0.5);
    double distance = z_height - z_max;
	double result = max_vmr*exp(-pow(distance, 2)/(2.0*pow(sigma, 2)));
	return result;
}

/*
Condensate properties
---------------------
*/

double c_v_cond(int solid_or_liquid, int subcategory, double temp)
{
	/*
	This function returns c_v of condensates.
	*/
	
	// solid_or_liquid == 0: ice
	// solid_or_liquid == 1: liquid water
	
    double result;
    if (solid_or_liquid == 0)
    {
        result = 2060.0;
    }
    if (solid_or_liquid == 1)
    {
        result = 4184.0;
    }
    return result;
}

double c_p_cond(int solid_or_liquid, int subcategory, double temp)
{
	/*
	This function returns c_p of condensates.
	*/
	
	// solid_or_liquid == 0: ice
	// solid_or_liquid == 1: liquid water
	
    double result;
    if (solid_or_liquid == 0)
    {
        result = 2060.0;
    }
    if (solid_or_liquid == 1)
    {
        result = 4184.0;
    }
    return result;
}

double phase_trans_heat(int direction, double temperature)
{
	/*
	This function calculates the phase transition heat.
	*/
	
    /*
    directions:
    0: gas to liquid
    1: gas to solid
    2: liquid to solid
    */
    
    double result = 0.0;
    if (direction == 0)
    {
    	result = enthalpy_evaporation(temperature);
	}
    if (direction == 1)
    {
        result = enthalpy_evaporation(temperature) + enthalpy_melting(temperature);
	}
    if (direction == 2)
    {
        result = enthalpy_melting(temperature);
	}
	
    return result;
}

double enthalpy_evaporation(double temperature)
{
	/*
	This function returns the enthalpy of evaporation depending on the temperature.
	It follows Pruppacher and Klett (2010), p. 97, Eq. (3-24a).
	*/
	
	// temperature in degrees Celsius
	double temp_c = temperature - T_0;
	
	// clipping values that are too extreme for these approximations
	if (temp_c < -20.0)
	{
		temp_c = -20.0;
	}
	if (temp_c > 40.0)
	{
		temp_c = 40.0;
	}
	
	double result = 597.3 - 0.561*temp_c;
	
	// unit conversion
	return 4186.8*result;
}

double enthalpy_melting(double temperature)
{
	/*
	This function returns the enthalpy of melting depending on the temperature.
	It follows Pruppacher and Klett (2010), p. 97, Eq. (3-26).
	*/
	
	// temperature in degrees Celsius
	double temp_c = temperature - T_0;
	
	// clipping values that are too extreme for this approximation
	if (temp_c < -44.0)
	{
		temp_c = -44.0;
	}
	if (temp_c > 0.0)
	{
		temp_c = 0.0;
	}
	
	double result = 79.7 - 0.12*temp_c - 8.0481e-2*pow(temp_c, 2) - 3.2376e-3*pow(temp_c, 3) - 4.2553e-5*pow(temp_c, 4);
	
	// unit conversion
	return 4186.8*result;
}

double saturation_pressure_over_water(double temperature)
{
	/*
	This function returns the saturation pressure in Pa over liquid water as a function of the temperature in K.
	It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
	*/
    
    // calculating the temperature in degrees Celsius
    double temp_c = temperature - T_0;
    
    // clipping too low values for stability reasons
    if (temp_c < -60.0)
    {
    	temp_c = -60.0;
    }
    
    double result;
    
    if (temp_c > 0.0)
    {
    	result = exp(34.494 - 4924.99/(temp_c + 237.1))/pow(temp_c + 105.0, 1.57);
    }
    // for super-cooled water we use the formula for ice
    else
    {
    	result = saturation_pressure_over_ice(temperature);
    }
    
    return result;
}

double saturation_pressure_over_ice(double temperature)
{
	/*
	This function returns the saturation pressure in Pa over ice as a function of the temperature in K.
	It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
	*/
    
    // calculating the temperature in degrees Celsius
    double temp_c = temperature - T_0;
    
    // clipping too low values for stability reasons
    if (temp_c < -60.0)
    {
    	temp_c = -60.0;
    }
    // at temperatures > 0 degrees C ice cannot exist in equilibrium which is why this is clipped
    if (temp_c > 0.0)
    {
    	temp_c = 0.0;
    }
    
    double result;
    
	result = exp(43.494 - 6545.8/(temp_c + 278.0))/pow(temp_c + 868.0, 2);
    
    return result;
}

double rel_humidity(double abs_humidity, double temperature)
{
	/*
	This function returns the relative humidity (NOT in percent) as a function of the absolute humidity in kg/m^3 and the temperature in K.
	*/
	
	double vapour_pressure = abs_humidity*specific_gas_constants(1)*temperature;
	double saturation_pressure;
	if (temperature > T_0)
	{
		saturation_pressure = saturation_pressure_over_water(temperature);
	}
	if (temperature <= T_0)
	{
		saturation_pressure = saturation_pressure_over_ice(temperature);
	}
	double result = vapour_pressure/saturation_pressure;
	return result;
}








