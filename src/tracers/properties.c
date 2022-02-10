/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains look-up functions for properties of the atmosphere.
*/

#include <math.h>

double c_v_cond(int solid_or_liquid, int subcategory, double temp)
{
	/*
	This function returns c_v of condensates.
	*/
    double result;
    if (solid_or_liquid == 0)
    {
        result = 2060;
    }
    if (solid_or_liquid == 1)
    {
        result = 4184;
    }
    return result;
}

double c_p_cond(int solid_or_liquid, int subcategory, double temp)
{
	/*
	This function returns c_p of condensates.
	*/
    double result;
    if (solid_or_liquid == 0)
    {
        result = 2060;
    }
    if (solid_or_liquid == 1)
    {
        result = 4184;
    }
    return result;
}


// gas quantities
// ------------------------

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

double mean_particle_masses_gas_lookup(int gas_constituent_id)
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

double spec_heat_capacities_v_gas_lookup(int gas_constituent_id)
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

double spec_heat_capacities_p_gas_lookup(int gas_constituent_id)
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

double specific_gas_constants_lookup(int gas_constituent_id)
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
	double sigma = fwhm/pow(8*log(2), 0.5);
    double distance = z_height - z_max;
	double result = max_vmr*exp(-pow(distance, 2)/(2*pow(sigma, 2)));
	return result;
}




