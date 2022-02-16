/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains functions calculating everything related to water vapour.
*/

#include <math.h>
#include "../game_constants.h"
#include "../game_types.h"
#include "constituents.h"

double saturation_pressure_over_water(double temperature)
{
	/*
	This function returns the saturation pressure in Pa over liquid water as a function of the temperature in K.
	*/
    double temp_c = temperature - T_0;
    double result = 100.0*6.112*exp(17.62*temp_c/(243.12 + temp_c));
    return result;
}

double saturation_pressure_over_ice(double temperature)
{
	/*
	This function returns the saturation pressure in Pa over ice as a function of the temperature in K.
	*/
    double temp_c = temperature - T_0;
    double result = 100.0*6.112*exp(22.46*temp_c/(272.62 + temp_c));
    return result;
}

double water_vapour_density_from_rel_humidity(double rel_humidity, double temperature, double density)
{
	/*
	This function returns the absolute humidity (water vapour density) in kg/m^3 as a function of the relative humidity (NOT in percent),
	the temperature in K and the density in kg/m^3.
	*/
    double water_vapour_density = rel_humidity*saturation_pressure_over_water(temperature)/(specific_gas_constants(1)*temperature);
    return water_vapour_density;
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






