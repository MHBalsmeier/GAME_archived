/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
The vertical advection of horizontal momentum is organized here.
*/

#include "../enum_and_typedefs.h"
#include <stdlib.h>
#include <stdio.h>
#include "../diagnostics/diagnostics.h"
#include <omp.h>
#include "../spatial_operators/spatial_operators.h"

int thomas_algorithm(double [], double [], double [], double [], double [], double [], double [], int);

int three_band_solver_hor(State *state_0, State *state_p1, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Implicit vertical advection of vertical momentum (Crank-Nicolson).
	Procedure derived in Kompendium.
	The algorithm follows https://de.wikipedia.org/wiki/Thomas-Algorithmus .
	*/
	double *a_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *b_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_prime_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_prime_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *vertical_velocity = malloc(NO_OF_LAYERS*sizeof(double));
	double *solution_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double delta_z, delta_v, dvdz;
	int i;
	for (i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		// diagnozing the vertical velocity
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			recov_hor_ver_pri(state_0 -> velocity_gas, j, i, &vertical_velocity[j], grid);
		}
		// filling up the original vectors
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			if (j == NO_OF_LAYERS - 2)
				delta_z = grid -> z_vector[NO_OF_VECTORS_V + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_VECTORS_V + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			else
				delta_z = grid -> z_vector[NO_OF_VECTORS_V + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_VECTORS_V + (j + 2)*NO_OF_VECTORS_PER_LAYER + i];
			a_vector[j] = vertical_velocity[j + 1]*delta_t/(2*delta_z);
			if (j == 0)
				delta_z = grid -> z_vector[NO_OF_VECTORS_V + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_VECTORS_V + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			else
				delta_z = grid -> z_vector[NO_OF_VECTORS_V + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_VECTORS_V + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			c_vector[j] = -vertical_velocity[j]*delta_t/(2*delta_z);
		}
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			if (j == 0)
			{
				delta_z = grid -> z_vector[NO_OF_VECTORS_V + j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_VECTORS_V + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
				b_vector[j] = 1 + delta_t*vertical_velocity[j]/(2*delta_z);
				delta_v = state_0 -> velocity_gas[NO_OF_VECTORS_V + j*NO_OF_VECTORS_PER_LAYER + i] - state_0 -> velocity_gas[NO_OF_VECTORS_V + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			}
			else if (j == NO_OF_LAYERS - 1)
			{
				delta_z = grid -> z_vector[NO_OF_VECTORS_V + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_VECTORS_V + j*NO_OF_VECTORS_PER_LAYER + i];
				b_vector[j] = 1 - delta_t*vertical_velocity[j]/(2*delta_z);
				delta_v = state_0 -> velocity_gas[NO_OF_VECTORS_V + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - state_0 -> velocity_gas[NO_OF_VECTORS_V + j*NO_OF_VECTORS_PER_LAYER + i];
			}
			else
			{
				b_vector[j] = 1;
				delta_z = grid -> z_vector[NO_OF_VECTORS_V + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_VECTORS_V + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
				delta_v = state_0 -> velocity_gas[NO_OF_VECTORS_V + (j - 1)*NO_OF_VECTORS_PER_LAYER + i] - state_0 -> velocity_gas[NO_OF_VECTORS_V + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			}
			dvdz = delta_v/delta_z;
			d_vector[j] = state_0 -> velocity_gas[NO_OF_VECTORS_V + j*NO_OF_VECTORS_PER_LAYER + i] + delta_t*state_tendency -> velocity_gas[NO_OF_VECTORS_V + j*NO_OF_VECTORS_PER_LAYER + i] - vertical_velocity[j]*delta_t/2*dvdz;
		}
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_p1 -> velocity_gas[NO_OF_VECTORS_V + j*NO_OF_VECTORS_PER_LAYER + i] = solution_vector[j];
		}
	}
	free(solution_vector);
	free(vertical_velocity);
	free(a_vector);
	free(b_vector);
	free(c_vector);
	free(d_vector);
	free(c_prime_vector);
	free(d_prime_vector);
	return 0;
}

int three_band_solver_ver_vel(State *state_0, State *state_p1, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Implicit vertical advection of vertical momentum (Euler).
	Procedure derived in Kompendium.
	The algorithm follows https://de.wikipedia.org/wiki/Thomas-Algorithmus .
	*/
	double *a_vector = malloc((NO_OF_LAYERS - 2)*sizeof(double));
	double *b_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *c_vector = malloc((NO_OF_LAYERS - 2)*sizeof(double));
	double *d_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *c_prime_vector = malloc((NO_OF_LAYERS - 2)*sizeof(double));
	double *d_prime_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *vertical_velocity = malloc(NO_OF_LAYERS*sizeof(double));
	double *solution_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double delta_z;
	int i;
	for (i = 0; i < NO_OF_VECTORS_V; ++i)
	{
		// extracting the vertical velocity
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			vertical_velocity[j] = state_0 -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER];
			vertical_velocity[j] = 0;
		}
		vertical_velocity[NO_OF_LAYERS - 1] = state_p1 -> velocity_gas[i + NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER];
		vertical_velocity[NO_OF_LAYERS - 1] = 0;
		// filling up the original vectors
		for (int j = 0; j < NO_OF_LAYERS - 2; ++j)
		{
			delta_z = grid -> z_vector[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 3)*NO_OF_VECTORS_PER_LAYER + i];
			a_vector[j] = vertical_velocity[j + 1]*delta_t/delta_z;
		}
		delta_z = grid -> z_vector[i] - grid -> z_vector[2*NO_OF_VECTORS_PER_LAYER + i];
		c_vector[0] = -vertical_velocity[0]*delta_t/delta_z;
		for (int j = 1; j < NO_OF_LAYERS - 2; ++j)
		{
			delta_z = grid -> z_vector[j*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[(j + 2)*NO_OF_VECTORS_PER_LAYER + i];
			c_vector[j] = -vertical_velocity[j]*delta_t/delta_z;
		}
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			b_vector[j] = 1;
			d_vector[j] = state_0 -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] + delta_t*state_tendency -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
		}
		delta_z = grid -> z_vector[(NO_OF_LAYERS - 2)*NO_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i];
		d_vector[NO_OF_LAYERS - 2] += delta_t*vertical_velocity[NO_OF_LAYERS - 2]*vertical_velocity[NO_OF_LAYERS - 1]/delta_z;
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS - 1);
		state_p1 -> velocity_gas[i] = 0;
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			state_p1 -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i] = solution_vector[j];
		}
	}
	free(solution_vector);
	free(vertical_velocity);
	free(a_vector);
	free(b_vector);
	free(c_vector);
	free(d_vector);
	free(c_prime_vector);
	free(d_prime_vector);
	return 0;
}

int three_band_solver_ver_den_dry(State *state_0, State *state_p1, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Implicit vertical advection of dry mass (Euler).
	Procedure derived in Kompendium.
	The algorithm follows https://de.wikipedia.org/wiki/Thomas-Algorithmus .
	*/
	double *a_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *b_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_prime_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_prime_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *vertical_flux_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *solution_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double area;
	int i;
	for (i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// diagnozing the vertical flux
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			if (j + 1 >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
				vertical_contravariant_normalized(state_p1 -> velocity_gas, j + 1, i, grid, &vertical_flux_vector[j]);
			else
				vertical_flux_vector[j] = state_p1 -> velocity_gas[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			vertical_flux_vector[j] = area*vertical_flux_vector[j];
		}
		// filling up the original vectors
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
		}
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			if (j == 0)
			{
				b_vector[j] = 1 - delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[0];
			}
			else if (j == NO_OF_LAYERS - 1)
			{
				b_vector[j] = 1 + delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j - 1];
			}
			else
			{
				b_vector[j] = 1 + delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*(vertical_flux_vector[j - 1] - vertical_flux_vector[j]);
			}
			d_vector[j] = state_0 -> density_dry[j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> density_dry[j*NO_OF_SCALARS_H + i];
		}
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_p1 -> density_dry[j*NO_OF_SCALARS_H + i] = solution_vector[j];
		}
	}
	free(solution_vector);
	free(vertical_flux_vector);
	free(a_vector);
	free(b_vector);
	free(c_vector);
	free(d_vector);
	free(c_prime_vector);
	free(d_prime_vector);
	return 0;
}

int three_band_solver_ver_entropy_gas(State *state_0, State *state_p1, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Implicit vertical advection of dry mass (Euler).
	Procedure derived in Kompendium.
	The algorithm follows https://de.wikipedia.org/wiki/Thomas-Algorithmus .
	*/
	double *a_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *b_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_prime_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_prime_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *vertical_flux_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *solution_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double area;
	int i;
	for (i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// diagnozing the vertical flux
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			if (j + 1 >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
				vertical_contravariant_normalized(state_p1 -> velocity_gas, j + 1, i, grid, &vertical_flux_vector[j]);
			else
				vertical_flux_vector[j] = state_p1 -> velocity_gas[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
			vertical_flux_vector[j] = area*vertical_flux_vector[j];
		}
		// filling up the original vectors
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			a_vector[j] = delta_t/(2*grid -> volume[i + (j + 1)*NO_OF_SCALARS_H])*vertical_flux_vector[j];
			c_vector[j] = -delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j];
		}
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			if (j == 0)
			{
				b_vector[j] = 1 - delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[0];
			}
			else if (j == NO_OF_LAYERS - 1)
			{
				b_vector[j] = 1 + delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*vertical_flux_vector[j - 1];
			}
			else
			{
				b_vector[j] = 1 + delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*(vertical_flux_vector[j - 1] - vertical_flux_vector[j]);
			}
			d_vector[j] = state_0 -> entropy_gas[j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> entropy_gas[j*NO_OF_SCALARS_H + i];
		}
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_p1 -> entropy_gas[j*NO_OF_SCALARS_H + i] = solution_vector[j];
		}
	}
	free(solution_vector);
	free(vertical_flux_vector);
	free(a_vector);
	free(b_vector);
	free(c_vector);
	free(d_vector);
	free(c_prime_vector);
	free(d_prime_vector);
	return 0;
}

int thomas_algorithm(double a_vector[], double b_vector[], double c_vector[], double d_vector[], double c_prime_vector[], double d_prime_vector[], double solution_vector[], int solution_length)
{
	// https://de.wikipedia.org/wiki/Thomas-Algorithmus
	c_prime_vector[0] = c_vector[0]/b_vector[0];
	for (int j = 1; j < solution_length - 1; ++j)
	{
		c_prime_vector[j] = c_vector[j]/(b_vector[j] - c_prime_vector[j - 1]*a_vector[j - 1]);
	}
	d_prime_vector[0] = d_vector[0]/b_vector[0];
	for (int j = 1; j < solution_length; ++j)
	{
		d_prime_vector[j] = (d_vector[j] - d_prime_vector[j - 1]*a_vector[j - 1])/(b_vector[j] - c_prime_vector[j - 1]*a_vector[j - 1]);
	}
	solution_vector[solution_length - 1] = d_prime_vector[solution_length - 1];
	for (int j = solution_length - 2; j >= 0; --j)
	{
		solution_vector[j] = d_prime_vector[j] - c_prime_vector[j]*solution_vector[j + 1];
	}
	return 0;
}





