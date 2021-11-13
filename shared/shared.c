/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

#include <math.h>
#include "geos95.h"
#define OMEGA (7.292115e-5)
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define P_0 100000.0

const double T_SFC = 273.15 + 15;
const double TEMP_GRADIENT = -0.65/100;
// constants that are specific to the ICAO standard atmosphere
const double P_0_STANDARD = 101325;
const double TROPO_HEIGHT_STANDARD = 11e3;
const double INVERSE_HEIGHT_STANDARD = 20e3;
const double TEMP_GRADIENT_INV_STANDARD = 0.1/100;

double standard_temp(double z_height)
{
    // temperature in the standard atmosphere
    const double TROPO_TEMP_STANDARD = T_SFC + TROPO_HEIGHT_STANDARD*TEMP_GRADIENT;
    double temperature;
    if (z_height < TROPO_HEIGHT_STANDARD)
    {
        temperature = T_SFC + z_height*TEMP_GRADIENT;
    }
    else if (z_height < INVERSE_HEIGHT_STANDARD)
    {
        temperature = TROPO_TEMP_STANDARD;
    }
    else
    {
    	temperature = TROPO_TEMP_STANDARD + TEMP_GRADIENT_INV_STANDARD*(z_height - INVERSE_HEIGHT_STANDARD);
    }
    return temperature;
}

double standard_pres(double z_height)
{
    // pressure in the standard atmosphere
	const double G = 9.80616;
    const double TROPO_TEMP_STANDARD = T_SFC + TROPO_HEIGHT_STANDARD*TEMP_GRADIENT;
    double pressure, pressure_at_inv_standard;
    if (z_height < TROPO_HEIGHT_STANDARD)
    {
        pressure = P_0_STANDARD*pow(1 + TEMP_GRADIENT*z_height/T_SFC, -G/(R_D*TEMP_GRADIENT));
    }
    else if (z_height < INVERSE_HEIGHT_STANDARD)
    {
        pressure = P_0_STANDARD*pow(1 + TEMP_GRADIENT*TROPO_HEIGHT_STANDARD/T_SFC, -G/(R_D*TEMP_GRADIENT))*exp(-G*(z_height - TROPO_HEIGHT_STANDARD)/(R_D*TROPO_TEMP_STANDARD));
    }
    else
    {
    	pressure_at_inv_standard = P_0_STANDARD*pow(1 + TEMP_GRADIENT*TROPO_HEIGHT_STANDARD/T_SFC, -G/(R_D*TEMP_GRADIENT))*exp(-G*(INVERSE_HEIGHT_STANDARD - TROPO_HEIGHT_STANDARD)/(R_D*TROPO_TEMP_STANDARD));
        pressure = pressure_at_inv_standard*pow(1 + TEMP_GRADIENT*(z_height - INVERSE_HEIGHT_STANDARD)/T_SFC, -G/(R_D*TEMP_GRADIENT));
    }
    return pressure;
}






















