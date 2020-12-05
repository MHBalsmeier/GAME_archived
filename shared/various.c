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

int lu_5band_solver(double a_vector[], double b_vector[], double c_vector[], double l_vector[], double u_vector[], double d_vector[], double solution_vector[], int solution_length)
{
	// Here the system of linear equations Ax = d is solved, using the LU decomposition.
	// We reformulate this problem to LUx = d with a lower triangular matrix L and an upper triangular matrix U.
	// Determining the vectors l_1_vector and l_2_vector, which make up the matrix L.
	double a_vector_mod[solution_length - 1];
	double b_vector_mod[solution_length];
	double c_vector_mod[solution_length - 1];
	for (int i = 0; i < solution_length; ++i)
	{
		b_vector_mod[i] = b_vector[i];
	}
	for (int i = 0; i < solution_length - 1; ++i)
	{
		a_vector_mod[i] = a_vector[i];
		c_vector_mod[i] = c_vector[i];
	}
	double l_1_vector[solution_length - 1];
	double l_2_vector[solution_length - 1];
	for (int i = 0; i < solution_length - 2; ++i)
	{
		l_1_vector[i] = a_vector_mod[i]/b_vector_mod[i];
		b_vector_mod[i + 1] += -c_vector_mod[i]*l_1_vector[i];
		c_vector_mod[i + 1] += -u_vector[i]*l_1_vector[i];
		l_2_vector[i] = l_vector[i]/b_vector_mod[i];
		a_vector_mod[i + 1] += -c_vector[i]*l_2_vector[i];
		b_vector_mod[i + 2] += -u_vector[i]*l_2_vector[i];
		if (i < solution_length - 3)
		{
			c_vector_mod[i + 2] += -u_vector[i + 1]*l_2_vector[i];
		}
	}
	l_1_vector[solution_length - 2] = a_vector_mod[solution_length - 2]/b_vector_mod[solution_length - 2];
	
	// Solving Ly = d.
	double y_vector[solution_length];
	y_vector[0] = d_vector[0];
	y_vector[1] = d_vector[1] - l_1_vector[0]*y_vector[0];
	for (int i = 2; i < solution_length; ++i)
	{
		y_vector[i] = d_vector[i] - l_1_vector[i - 1]*y_vector[i - 1] - l_2_vector[i - 2]*y_vector[i - 2];
	}
	
	// It is L^{-1}A = R:
	double r_0_vector[solution_length];
	double r_1_vector[solution_length - 1];
	r_0_vector[0] = b_vector[0];
	r_1_vector[0] = c_vector[0];
	r_0_vector[1] = -l_1_vector[0]*c_vector[0] + b_vector[1];
	for (int i = 2; i < solution_length; ++i)
	{
		r_0_vector[i] = -l_2_vector[i - 2]*u_vector[i - 2] - l_1_vector[i - 1]*c_vector[i - 1] + b_vector[i];
	}
	r_1_vector[1] = -l_1_vector[0]*u_vector[0] + c_vector[1];
	for (int i = 2; i < solution_length - 1; ++i)
	{
		r_1_vector[i] = -l_1_vector[i - 1]*u_vector[i - 1] + c_vector[i];
	}
	
	// Solving Rx = y.
	solution_vector[solution_length - 1] = y_vector[solution_length - 1]/r_0_vector[solution_length - 1];
	solution_vector[solution_length - 2] = (y_vector[solution_length - 2] - r_1_vector[solution_length - 2]*solution_vector[solution_length - 1])/r_0_vector[solution_length - 2];
	for (int i = solution_length - 3; i >= 0; --i)
	{
		solution_vector[i] = (y_vector[i] - r_1_vector[i]*solution_vector[i + 1] - u_vector[i]*solution_vector[i + 2])/r_0_vector[i];
	}
	return 0;
}



