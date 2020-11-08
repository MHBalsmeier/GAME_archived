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
    {
        phi_bg = T_0*G/GAMMA*(1 - pow(eta, R_D*GAMMA/G));
    }
    else
    {
        phi_bg = T_0*G/GAMMA*(1 - pow(eta, R_D*GAMMA/G)) - R_D*DELTA_T*((log(eta/ETA_T) + 137.0/60)*pow(ETA_T, 5) - 5*eta*pow(ETA_T, 4) + 5*pow(ETA_T, 3)*pow(eta, 2) - 10.0/3*pow(ETA_T, 2)*pow(eta, 3) + 5.0/4*ETA_T*pow(eta, 4) - 1.0/5*pow(eta, 5));
    }
    phi_perturb = U_0*pow(cos(eta_v), 1.5)*((-2*pow(sin(lat), 6)*(pow(cos(lat), 2) + 1.0/3) + 10.0/63)*U_0*pow(cos(eta_v), 1.5) + RADIUS*OMEGA*(8.0/5*pow(cos(lat), 3)*(pow(sin(lat), 2) + 2.0/3) - M_PI/4));
    phi = phi_bg + phi_perturb;
    z = phi/G;
    *result = z;
    return 0;
}
