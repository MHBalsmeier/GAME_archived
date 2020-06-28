/*
This source file is part of the Global Geophysical Modeling Frame (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
The vertical advection of horizontal momentum is organized here.
*/

#include "../enum_and_typedefs.h"
#include <stdlib.h>
#include <stdio.h>
#include "../diagnostics/diagnostics.h"

int three_band_solver(State *state_0, State *state_p1, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	The algorithm follows https://de.wikipedia.org/wiki/Thomas-Algorithmus and Kompendium. Notation follows the Wikipedia article.
	*/
	double *a_vector = malloc((NUMBER_OF_LAYERS - 1)*sizeof(double));
	double *b_vector = malloc(NUMBER_OF_LAYERS*sizeof(double));
	double *c_vector = malloc((NUMBER_OF_LAYERS - 1)*sizeof(double));
	double *d_vector = malloc(NUMBER_OF_LAYERS*sizeof(double));
	double *c_prime_vector = malloc((NUMBER_OF_LAYERS - 1)*sizeof(double));
	double *d_prime_vector = malloc(NUMBER_OF_LAYERS*sizeof(double));
	double *vertical_velocity = malloc(NUMBER_OF_LAYERS*sizeof(double));
	double delta_z, delta_v, dvdz;
	for (int i = 0; i < NUMBER_OF_VECTORS_H; ++i)
	{
		// diagnozing the vertical velocity
		for (int j = 0; j < NUMBER_OF_LAYERS; ++j)
		{
			recov_hor_ver_pri(state_0 -> velocity_gas, j, i, &vertical_velocity[j], grid);
		}
		// filling up the original vectors
		for (int j = 0; j < NUMBER_OF_LAYERS - 1; ++j)
		{
			if (j == NUMBER_OF_LAYERS - 2)
				delta_z = grid -> z_vector[NUMBER_OF_VECTORS_V + j*NUMBER_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NUMBER_OF_VECTORS_V + (j + 1)*NUMBER_OF_VECTORS_PER_LAYER + i];
			else
				delta_z = grid -> z_vector[NUMBER_OF_VECTORS_V + j*NUMBER_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NUMBER_OF_VECTORS_V + (j + 2)*NUMBER_OF_VECTORS_PER_LAYER + i];
			a_vector[j] = -vertical_velocity[j]*delta_t/(2*delta_z);
			if (j == 0)
				delta_z = grid -> z_vector[NUMBER_OF_VECTORS_V + j*NUMBER_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NUMBER_OF_VECTORS_V + (j + 1)*NUMBER_OF_VECTORS_PER_LAYER + i];
			else
				delta_z = grid -> z_vector[NUMBER_OF_VECTORS_V + (j - 1)*NUMBER_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NUMBER_OF_VECTORS_V + (j + 1)*NUMBER_OF_VECTORS_PER_LAYER + i];
			c_vector[j] = vertical_velocity[j]*delta_t/(2*delta_z);
		}
		for (int j = 0; j < NUMBER_OF_LAYERS; ++j)
		{
			if (j == 0)
			{
				delta_z = grid -> z_vector[NUMBER_OF_VECTORS_V + j*NUMBER_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NUMBER_OF_VECTORS_V + (j + 1)*NUMBER_OF_VECTORS_PER_LAYER + i];
				b_vector[j] = 1 - delta_t*vertical_velocity[j]/(2*delta_z);
				delta_v = state_0 -> velocity_gas[NUMBER_OF_VECTORS_V + j*NUMBER_OF_VECTORS_PER_LAYER + i] - state_0 -> velocity_gas[NUMBER_OF_VECTORS_V + (j + 1)*NUMBER_OF_VECTORS_PER_LAYER + i];
			}
			else if (j == NUMBER_OF_LAYERS - 1)
			{
				delta_z = grid -> z_vector[NUMBER_OF_VECTORS_V + (j - 1)*NUMBER_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NUMBER_OF_VECTORS_V + j*NUMBER_OF_VECTORS_PER_LAYER + i];
				b_vector[j] = 1 + delta_t*vertical_velocity[j]/(2*delta_z);
				delta_v = state_0 -> velocity_gas[NUMBER_OF_VECTORS_V + (j - 1)*NUMBER_OF_VECTORS_PER_LAYER + i] - state_0 -> velocity_gas[NUMBER_OF_VECTORS_V + j*NUMBER_OF_VECTORS_PER_LAYER + i];
			}
			else
			{
				b_vector[j] = 1;
				delta_z = grid -> z_vector[NUMBER_OF_VECTORS_V + (j - 1)*NUMBER_OF_VECTORS_PER_LAYER + i] - grid -> z_vector[NUMBER_OF_VECTORS_V + (j + 1)*NUMBER_OF_VECTORS_PER_LAYER + i];
				delta_v = state_0 -> velocity_gas[NUMBER_OF_VECTORS_V + (j - 1)*NUMBER_OF_VECTORS_PER_LAYER + i] - state_0 -> velocity_gas[NUMBER_OF_VECTORS_V + (j + 1)*NUMBER_OF_VECTORS_PER_LAYER + i];
			}
			dvdz = delta_v/delta_z;
			d_vector[j] = state_0 -> velocity_gas[NUMBER_OF_VECTORS_V + j*NUMBER_OF_VECTORS_PER_LAYER + i] + delta_t*state_tendency -> velocity_gas[NUMBER_OF_VECTORS_V + j*NUMBER_OF_VECTORS_PER_LAYER + i] - vertical_velocity[j]*delta_t/2*dvdz;
		}
		// modified vectors
		c_prime_vector[0] = c_vector[0]/b_vector[0];
		for (int j = 1; j < NUMBER_OF_LAYERS - 1; ++j)
		{
			c_prime_vector[j] = c_vector[j]/(b_vector[j] - c_prime_vector[j - 1]*a_vector[j]);
		}
		d_prime_vector[0] = d_vector[0]/b_vector[0];
		for (int j = 1; j < NUMBER_OF_LAYERS; ++j)
		{
			d_prime_vector[j] = (d_vector[j] - d_prime_vector[j - 1]*a_vector[j])/(b_vector[j] - c_prime_vector[j - 1]*a_vector[j]);
		}
		// finally, the solution is done here
		state_p1 -> velocity_gas[NUMBER_OF_VECTORS_V + (NUMBER_OF_LAYERS - 1)*NUMBER_OF_VECTORS_PER_LAYER + i] = d_prime_vector[NUMBER_OF_LAYERS - 1];
		for (int j = NUMBER_OF_LAYERS - 2; j >= 0; --j)
		{
			state_p1 -> velocity_gas[NUMBER_OF_VECTORS_V + j*NUMBER_OF_VECTORS_PER_LAYER + i] = d_prime_vector[j] - c_prime_vector[j]*state_p1 -> velocity_gas[NUMBER_OF_VECTORS_V + (j + 1)*NUMBER_OF_VECTORS_PER_LAYER + i];
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












