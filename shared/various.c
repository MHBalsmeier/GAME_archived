/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
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

int find_z_from_p(double lat, double p, double *result)
{
	const double ETA_0 = 0.252;
	const double ETA_T = 0.2;
	const double DELTA_T = 4.8e5;
	const double U_0 = 35;
	double GAMMA = 0.005;
	const double G = 9.80616;
	double T_0 = 288;
	// this function converts a preessure value into a geomtrical height (as a function of latitude) for the JW test
    double z;
    double eta = p/P_0;
    double phi, phi_bg, phi_perturb;
    double eta_v = (eta - ETA_0)*M_PI/2;
    if (eta >= ETA_T)
        phi_bg = T_0*G/GAMMA*(1 - pow(eta, R_D*GAMMA/G));
    else
        phi_bg = T_0*G/GAMMA*(1 - pow(eta, R_D*GAMMA/G)) - R_D*DELTA_T*((log(eta/ETA_T) + 137.0/60)*pow(ETA_T, 5) - 5*eta*pow(ETA_T, 4) + 5*pow(ETA_T, 3)*pow(eta, 2) - 10.0/3*pow(ETA_T, 2)*pow(eta, 3) + 5.0/4*ETA_T*pow(eta, 4) - 1.0/5*pow(eta, 5));
    phi_perturb = U_0*pow(cos(eta_v), 1.5)*((-2*pow(sin(lat), 6)*(pow(cos(lat), 2) + 1.0/3) + 10.0/63)*U_0*pow(cos(eta_v), 1.5) + RADIUS*OMEGA*(8.0/5*pow(cos(lat), 3)*(pow(sin(lat), 2) + 2.0/3) - M_PI/4));
    phi = phi_bg + phi_perturb;
    z = phi/G;
    *result = z;
    return 0;
}

// thermodynamic quantities
// ------------------------

/*
gaseous constituents IDs:
0: dry air
1: water vapour
*/

double entropy_constants_gas_lookup(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 2429487178047751925300627872548148580712448.000000; // ((MEAN_MASS_D*exp(5.0/3))/(3*M_PI*H_BAR*H_BAR))
	}
	if (gas_constituent_id == 1)
	{
		result = 1511084890012154487904341578321985168998400.000000; // ((MEAN_MASS_V*exp(5.0/3))/(3*M_PI*H_BAR*H_BAR))
	}
	return result;
}

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

