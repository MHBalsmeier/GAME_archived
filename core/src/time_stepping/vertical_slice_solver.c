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
		// modified vectors
		c_prime_vector[0] = c_vector[0]/b_vector[0];
		for (int j = 1; j < NO_OF_LAYERS - 1; ++j)
		{
			c_prime_vector[j] = c_vector[j]/(b_vector[j] - c_prime_vector[j - 1]*a_vector[j]);
		}
		d_prime_vector[0] = d_vector[0]/b_vector[0];
		for (int j = 1; j < NO_OF_LAYERS; ++j)
		{
			d_prime_vector[j] = (d_vector[j] - d_prime_vector[j - 1]*a_vector[j])/(b_vector[j] - c_prime_vector[j - 1]*a_vector[j]);
		}
		// finally, the solution is done here
		state_p1 -> velocity_gas[NO_OF_VECTORS_V + (NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + i] = d_prime_vector[NO_OF_LAYERS - 1];
		for (int j = NO_OF_LAYERS - 2; j >= 0; --j)
		{
			state_p1 -> velocity_gas[NO_OF_VECTORS_V + j*NO_OF_VECTORS_PER_LAYER + i] = d_prime_vector[j] - c_prime_vector[j]*state_p1 -> velocity_gas[NO_OF_VECTORS_V + (j + 1)*NO_OF_VECTORS_PER_LAYER + i];
			
		}
	}
	free(vertical_velocity);
	free(a_vector);
	free(b_vector);
	free(c_vector);
	free(d_vector);
	free(c_prime_vector);
	free(d_prime_vector);
	return 0;
}

int three_band_solver_ver(State *state_0, State *state_p1, State *state_tendency, double delta_t, Grid *grid)
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
		// modified vectors
		c_prime_vector[0] = c_vector[0]/b_vector[0];
		for (int j = 1; j < NO_OF_LAYERS - 2; ++j)
		{
			c_prime_vector[j] = c_vector[j]/(b_vector[j] - c_prime_vector[j - 1]*a_vector[j]);
		}
		d_prime_vector[0] = d_vector[0]/b_vector[0];
		for (int j = 1; j < NO_OF_LAYERS - 1; ++j)
		{
			d_prime_vector[j] = (d_vector[j] - d_prime_vector[j - 1]*a_vector[j])/(b_vector[j] - c_prime_vector[j - 1]*a_vector[j]);
		}
		// finally, the solution is done here
		state_p1 -> velocity_gas[(NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + i] = d_prime_vector[NO_OF_LAYERS - 2];
		for (int j = NO_OF_LAYERS - 2; j >= 1; --j)
		{
			state_p1 -> velocity_gas[j*NO_OF_VECTORS_PER_LAYER + i] = d_prime_vector[j - 1] - c_prime_vector[j]*state_p1 -> velocity_gas[(j + 1)*NO_OF_VECTORS_PER_LAYER + i];
		}
		state_p1 -> velocity_gas[i] = 0;
	}
	free(vertical_velocity);
	free(a_vector);
	free(b_vector);
	free(c_vector);
	free(d_vector);
	free(c_prime_vector);
	free(d_prime_vector);
	return 0;
}










