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






















