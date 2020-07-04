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
#include "atmostracers.h"

int thomas_algorithm(double [], double [], double [], double [], double [], double [], double [], int);
int candiate_has_converged(State *, int, double []);
int convergence_checker(double [], double []);
int modify_candidate(double [], double [], double [], double);
int calc_implicit_tendencies(double [], double [], Grid *, int, double []);

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
			recov_hor_ver_pri(state_p1 -> velocity_gas, j, i, &vertical_velocity[j], grid);
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
	Can be considered a preconditioner for the vertical sound wave solver.
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
			vertical_velocity[j] = state_p1 -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER];
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
	Can be considered a preconditioner for the vertical sound wave solver.
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
	Implicit vertical advection of the entropy of the gas phase (Euler).
	Procedure derived in Kompendium.
	The algorithm follows https://de.wikipedia.org/wiki/Thomas-Algorithmus .
	Can be considered a preconditioner for the vertical sound wave solver.
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

int three_band_solver_ver_tracers(State *state_0, State *state_p1, State *state_tendency, double delta_t, Grid *grid)
{
	/*
	Implicit vertical advection of tracers (Euler).
	Procedure derived in Kompendium.
	The algorithm follows https://de.wikipedia.org/wiki/Thomas-Algorithmus .
	*/
	double *a_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *b_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *c_prime_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *d_prime_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *solution_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *vertical_flux_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double area;
	int i;
	for (i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		for (int k = 0; k < NO_OF_CONDENSATED_TRACERS; ++k)
		{
			// tracer masses
			// diagnozing the vertical flux
			for (int j = 0; j < NO_OF_LAYERS; ++j)
			{
				area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
				if (j + 1 >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
				{
					vertical_contravariant_normalized(state_p1 -> velocity_gas, j + 1, i, grid, &vertical_flux_vector[j]);
					// small metric terms neglected here
					vertical_flux_vector[j] -= ret_sink_velocity(i, 0, 0.001);
				}
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
				else
				{
					b_vector[j] = 1 + delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*(vertical_flux_vector[j - 1] - vertical_flux_vector[j]);
				}
				d_vector[j] = state_0 -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
			}
			thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
			for (int j = 0; j < NO_OF_LAYERS; ++j)
			{
				state_p1 -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = solution_vector[j];
				// limiter
				if (state_p1 -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] < 0)
					state_p1 -> tracer_densities[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = 0;
			}
			// tracer temberature densities
			// diagnozing the vertical flux
			for (int j = 0; j < NO_OF_LAYERS; ++j)
			{
				area = grid -> area[NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
				if (j + 1 >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
				{
					vertical_contravariant_normalized(state_p1 -> velocity_gas, j + 1, i, grid, &vertical_flux_vector[j]);
					// small metric term negelcted here
					vertical_flux_vector[j] -= ret_sink_velocity(i, 0, 0.001);
				}
				else
					vertical_flux_vector[j] = state_p1 -> velocity_gas[k*NO_OF_SCALARS + NO_OF_VECTORS_PER_LAYER + j*NO_OF_VECTORS_PER_LAYER + i];
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
				else
				{
					b_vector[j] = 1 + delta_t/(2*grid -> volume[i + j*NO_OF_SCALARS_H])*(vertical_flux_vector[j - 1] - vertical_flux_vector[j]);
				}
				d_vector[j] = state_0 -> tracer_density_temperatures[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> tracer_density_temperatures[k*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
			}
			thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
			for (int j = 0; j < NO_OF_LAYERS; ++j)
			{
				state_p1 -> tracer_density_temperatures[j*NO_OF_SCALARS_H + i] = solution_vector[j];
				// limiter
				if (state_p1 -> tracer_density_temperatures[j*NO_OF_SCALARS_H + i] < 0)
					state_p1 -> tracer_density_temperatures[j*NO_OF_SCALARS_H + i] = 0;
			}
		}
		// water vaour
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
			d_vector[j] = state_0 -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] + delta_t*state_tendency -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i];
		}
		thomas_algorithm(a_vector, b_vector, c_vector, d_vector, c_prime_vector, d_prime_vector, solution_vector, NO_OF_LAYERS);
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			state_p1 -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = solution_vector[j];
			if (state_p1 -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] < 0)
				state_p1 -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + j*NO_OF_SCALARS_H + i] = 0;
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

int vertical_sound_wave_solver(State *state_0, State *state_p1, State *state_tendency, double delta_t, Grid *grid, Vector_field pressure_gradient_acc, int step_counter)
{
	/*
	This is the iterative solver required for stabilizing vertically propagating sound waves.
	*/
	double *candidate_vector = malloc((2*NO_OF_LAYERS + NO_OF_LAYERS - 1)*sizeof(double));
	double *candidate_vector_result = malloc((2*NO_OF_LAYERS + NO_OF_LAYERS - 1)*sizeof(double));
	double *implicit_tendency = malloc((2*NO_OF_LAYERS + NO_OF_LAYERS - 1)*sizeof(double));
	int has_converged, number_of_iterations;
	double *vertical_velocity_gradient = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{	
		for (int j = 0; j < NO_OF_LAYERS; ++j)
		{
			candidate_vector[3*j] = state_p1 -> density_dry[i + j*NO_OF_SCALARS_H];
			candidate_vector[3*j + 1] = state_p1 -> entropy_gas[i + j*NO_OF_SCALARS_H];
		}
		for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
		{
			candidate_vector[3*j + 2] = state_p1 -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER];
		}
		has_converged = 0;
		number_of_iterations = 0;
		while (has_converged == 0)
		{
			calc_implicit_tendencies(candidate_vector, implicit_tendency, grid, i, vertical_velocity_gradient);
			for (int j = 0; j < NO_OF_LAYERS; ++j)
			{
				candidate_vector_result[3*j] = state_0 -> density_dry[i + j*NO_OF_SCALARS_H] + delta_t*state_tendency -> density_dry[i + j*NO_OF_SCALARS_H] + delta_t*implicit_tendency[3*j];
				candidate_vector_result[3*j + 1] = state_0 -> entropy_gas[i + j*NO_OF_SCALARS_H] + delta_t*state_tendency -> entropy_gas[i + j*NO_OF_SCALARS_H] + delta_t*implicit_tendency[3*j + 1];
				/*if (step_counter == 2)
				{
					printf("%d\n", number_of_iterations);
					printf("%lf\n", implicit_tendency[3*j]);
				}*/
			}
			for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
			{
				candidate_vector_result[3*j + 2] = state_0 -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER] + delta_t*(state_tendency -> velocity_gas[i + (j + 1)*NO_OF_VECTORS_PER_LAYER] - pressure_gradient_acc[i + (j + 1)*NO_OF_VECTORS_PER_LAYER]) + delta_t*implicit_tendency[3*j + 2];
			}
			has_converged = convergence_checker(candidate_vector_result, candidate_vector);
			if (has_converged == 0)
			{
				modify_candidate(candidate_vector_result, candidate_vector, vertical_velocity_gradient, delta_t);
			}
			++number_of_iterations;
		}
		candiate_has_converged(state_p1, i, candidate_vector);
	}
	free(vertical_velocity_gradient);
	free(implicit_tendency);
	free(candidate_vector);
	free(candidate_vector_result);
	return 0;
}

int calc_implicit_tendencies(double candidate_vector[], double implicit_tendency[], Grid *grid, int h_index, double vertical_velocity_gradient[])
{
	double *dry_mass_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *entropy_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double *vertical_velocity_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	for (int j = 0; j < NO_OF_LAYERS; ++j)
	{
		dry_mass_vector[j] = candidate_vector[3*j];
		entropy_vector[j] = candidate_vector[3*j + 1];
	}
	for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
	{
		vertical_velocity_vector[j] = candidate_vector[3*j + 2];
	}
	double *dry_mass_flux_density = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	scalar_times_vector_v_column(dry_mass_vector, vertical_velocity_vector, dry_mass_flux_density);
	double *entropy_flux_density = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	scalar_times_vector_v_column(entropy_vector, vertical_velocity_vector, entropy_flux_density);
	double *dry_mass_flux_density_divergence = malloc(NO_OF_LAYERS*sizeof(double));
	divv_v_columns(dry_mass_flux_density, dry_mass_flux_density_divergence, h_index, grid);
	double *entropy_flux_density_divergence = malloc(NO_OF_LAYERS*sizeof(double));
	divv_v_columns(entropy_flux_density, entropy_flux_density_divergence, h_index, grid);
	for (int j = 0; j < NO_OF_LAYERS; ++j)
	{
		implicit_tendency[3*j] = -dry_mass_flux_density_divergence[j];
		implicit_tendency[3*j + 1] = -entropy_flux_density_divergence[j];
	}
	double *temperature_vector = malloc(NO_OF_LAYERS*sizeof(double));
	double pot_temp_value, exner_pressure_value;
	for (int j = 0; j < NO_OF_LAYERS; ++j)
	{
    	pot_temp_value = pot_temp_diagnostics_single_value(entropy_vector[j], dry_mass_vector[j], 0.0, 0.0);
        exner_pressure_value = exner_pressure_diagnostics_single_value(dry_mass_vector[j], 0.0, pot_temp_value);
		temperature_vector[j] = temperature_diagnostics_single_value(exner_pressure_value, pot_temp_value);
	}
	double *pressure_gradient_acc_vector_0 = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	grad_v_scalar_column(temperature_vector, pressure_gradient_acc_vector_0, h_index, grid);
	double *pressure_gradient_acc_vector_1 = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	double *entropy_specific_vector = malloc(NO_OF_LAYERS*sizeof(double));
	for (int j = 0; j < NO_OF_LAYERS; ++j)
		entropy_specific_vector[j] = entropy_vector[j]/dry_mass_vector[j];
	grad_v_scalar_column(entropy_specific_vector, pressure_gradient_acc_vector_1, h_index, grid);
	free(entropy_specific_vector);
	scalar_times_vector_v_column(temperature_vector, pressure_gradient_acc_vector_1, pressure_gradient_acc_vector_1);
	double *pressure_gradient_acc_vector = malloc((NO_OF_LAYERS - 1)*sizeof(double));
	for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
	{
		pressure_gradient_acc_vector[j] = -C_D_P*pressure_gradient_acc_vector_0[j] + pressure_gradient_acc_vector_1[j];
	}
	grad_v_vector_column(vertical_velocity_vector, vertical_velocity_gradient, h_index, grid);
	for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
	{
		implicit_tendency[3*j + 2] = -vertical_velocity_vector[j]*vertical_velocity_gradient[j] + pressure_gradient_acc_vector[j];
	}
	free(pressure_gradient_acc_vector);
	free(pressure_gradient_acc_vector_0);
	free(pressure_gradient_acc_vector_1);
	free(dry_mass_vector);
	free(entropy_vector);
	free(vertical_velocity_vector);
	free(entropy_flux_density);
	free(dry_mass_flux_density);
	free(entropy_flux_density_divergence);
	free(dry_mass_flux_density_divergence);
	return 0;
}

int modify_candidate(double candidate_vector_result[], double candidate_vector[], double dwdz_vector[], double delta_t)
{
	double dwdz;
	for (int i = 0; i < NO_OF_LAYERS; ++i)
	{
		if (i == 0)
			dwdz = 0.5*dwdz_vector[0];
		else if (i == NO_OF_LAYERS - 1)
			dwdz = 0.5*dwdz_vector[i - 1];
		else
			dwdz = 0.5*(dwdz_vector[i - 1] + dwdz_vector[i]);
		candidate_vector[3*i + 0] = candidate_vector[3*i + 0] + (candidate_vector_result[3*i + 0] - candidate_vector[3*i + 0] - 0.5*delta_t*dwdz);
		candidate_vector[3*i + 1] = candidate_vector[3*i + 1] + (candidate_vector_result[3*i + 1] - candidate_vector[3*i + 1] - 0.5*delta_t*dwdz);
	}
	for (int i = 0; i < NO_OF_LAYERS - 1; ++i)
	{
		candidate_vector[3*i + 2] = candidate_vector[3*i + 2] + (candidate_vector_result[3*i + 2] - candidate_vector[3*i + 2] - 0.5*delta_t*dwdz_vector[i]);
	}
	return 0;
}

int convergence_checker(double candidate_vector_result[], double candidate_vector[])
{
	// this functions checks wether the vertical sound wave solver has converged
	int result = 1;
	double density_criterion = 0.001;
	double entropy_criterion = 0.1;
	double vertical_velocity_criterion = 0.001;
	for (int j = 0; j < NO_OF_LAYERS; ++j)
	{
		if (fabs(candidate_vector_result[3*j + 0] - candidate_vector[3*j + 0]) >= density_criterion)
			result = 0;
		if (fabs(candidate_vector_result[3*j + 1] - candidate_vector[3*j + 1]) >= entropy_criterion)
			result = 0;
	}
	for (int j = 0; j < NO_OF_LAYERS - 1; ++j)
	{
		if (fabs(candidate_vector_result[3*j + 2] - candidate_vector[3*j + 2]) >= vertical_velocity_criterion)
			result = 0;
	}
	return result;
}

int candiate_has_converged(State *state_p1, int h_index, double candidate_vector[])
{
	for (int i = 0; i < NO_OF_LAYERS; ++i)
	{
		state_p1 -> density_dry[h_index + i*NO_OF_SCALARS_H] = candidate_vector[3*i];
		state_p1 -> entropy_gas[h_index + i*NO_OF_SCALARS_H] = candidate_vector[3*i + 1];
	}
	for (int i = 0; i < NO_OF_LAYERS - 1; ++i)
	{
		state_p1 -> velocity_gas[h_index + (i + 1)*NO_OF_VECTORS_PER_LAYER] = candidate_vector[3*i + 2];
	}
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





