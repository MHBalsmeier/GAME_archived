/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, algebraic conversions and calculations of thermodynamic quantities of a moist atmosphere are collected.
indices as usual:
d:	dry
v:	water vapour
h:	humid
*/

#include <stdio.h>
#include <stdlib.h>
#include <atmostracers.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../spatial_operators/spatial_operators.h"
#include "thermodynamics.h"

int temperature_diagnostics(State *state, Grid *grid, Diagnostics *diagnostics)
{
	/*
	This function diagnoses the temperature of the gas phase.
	*/
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> temperature_gas[i] = (grid -> theta_bg[i] + state -> theta_pert[i])*(grid -> exner_bg[i] + state -> exner_pert[i]);
	}
	return 0;
}

double spec_heat_cap_diagnostics_v(State *state, int grid_point_index, Config *config)
{
	double rho_g = 0;
	int no_of_relevant_constituents = 0;
	if (config -> assume_lte == 0)
	{
		no_of_relevant_constituents = NO_OF_GASEOUS_CONSTITUENTS;
		rho_g = density_gas(state, grid_point_index);
	}
	if (config -> assume_lte == 1)
	{
		no_of_relevant_constituents = 1;
		rho_g = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + grid_point_index];
	}
	double result = 0;
	for (int i = 0; i < no_of_relevant_constituents; ++i)
	{
		result += state -> rho[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index]/rho_g*spec_heat_capacities_v_gas(i);
	}
	return result;
}

double spec_heat_cap_diagnostics_p(State *state, int grid_point_index, Config *config)
{
	double rho_g = 0;
	int no_of_relevant_constituents = 0;
	if (config -> assume_lte == 0)
	{
		no_of_relevant_constituents = NO_OF_GASEOUS_CONSTITUENTS;
		rho_g = density_gas(state, grid_point_index);
	}
	if (config -> assume_lte == 1)
	{
		no_of_relevant_constituents = 1;
		rho_g = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + grid_point_index];
	}
	double result = 0;
	for (int i = 0; i < no_of_relevant_constituents; ++i)
	{
		result += state -> rho[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index]/rho_g*spec_heat_capacities_p_gas(i);
	}
	return result;
}

double gas_constant_diagnostics(State *state, int grid_point_index, Config *config)
{
	double rho_g = 0;
	int no_of_relevant_constituents = 0;
	if (config -> assume_lte == 0)
	{
		no_of_relevant_constituents = NO_OF_GASEOUS_CONSTITUENTS;
		rho_g = density_gas(state, grid_point_index);
	}
	if (config -> assume_lte == 1)
	{
		no_of_relevant_constituents = 1;
		rho_g = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + grid_point_index];
	}
	double result = 0;
	for (int i = 0; i < no_of_relevant_constituents; ++i)
	{
		result += state -> rho[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index]/rho_g*specific_gas_constants(i);
	}
	return result;
}

double density_total(State *state, int grid_point_index)
{
	double result = 0;
	for (int i = 0; i < NO_OF_CONSTITUENTS; ++i)
	{
		result += state -> rho[i*NO_OF_SCALARS + grid_point_index];
	}
	return result;
}

double density_gas(State *state, int grid_point_index)
{
	double result = 0;
	for (int i = 0; i < NO_OF_GASEOUS_CONSTITUENTS; ++i)
	{
		result += state -> rho[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index];
	}
	return result;
}

double calc_micro_density(double density_macro, double condensates_density_sum)
{
	/*
	In a moist atmosphere one needs to distinguish between the densities with respect to the whole volume and the densities with respect to exclusively the gas phase.
	*/
	double result = density_macro/(1 - condensates_density_sum/RHO_WATER);
	if (result < 0)
	{
		printf("Error: microscopic density is negative.\n");
		printf("Aborting.\n");
		exit(1);
	}
	if (isnan(result))
	{
		printf("Error: microscopic density is nan.\n");
		printf("Aborting.\n");
		exit(1);
	}
	return result;
}

double calc_condensates_density_sum(int scalar_gridpoint_index, Mass_densities mass_densities)
{
	/*
	This is only needed for calculating the "micro densities".
	*/
	double result = 0;
	for (int i = 0; i < NO_OF_CONDENSED_CONSTITUENTS; ++i)
	{
		result += mass_densities[i*NO_OF_SCALARS + scalar_gridpoint_index];
	}
	if (result < 0)
	{
		printf("Error: condensates_density_sum is negative.\n");
		printf("Aborting.\n");
		exit(1);
	}
	if (result >= RHO_WATER)
	{
		printf("Error: condensates_density_sum >= RHO_WATER.\n");
		printf("Aborting.\n");
		exit(1);
	}
	return result;
}

double calc_diffusion_coeff(double temperature, double density)
{
	/*
	This function calculates the molecular diffusion coefficient according to the kinetic gas theory.
	*/
	// these things are hardly ever modified
	double particle_radius = 130e-12;
	double particle_mass = mean_particle_masses_gas(0);
	// actual calculation
    double thermal_velocity = sqrt(8*K_B*temperature/(M_PI*particle_mass));
    double particle_density = density/particle_mass;
    double cross_section = 4*M_PI*pow(particle_radius, 2);
    double mean_free_path = 1/(sqrt(2)*particle_density*cross_section);
    double result = 1.0/3*thermal_velocity*mean_free_path;
    return result;
}


double mean_particle_masses_gas(int gas_constituent_id)
{
	// binding to atmostracers
	return mean_particle_masses_gas_lookup(gas_constituent_id);
}

double spec_heat_capacities_v_gas(int gas_constituent_id)
{
	// binding to atmostracers
	return spec_heat_capacities_v_gas_lookup(gas_constituent_id);
}

double spec_heat_capacities_p_gas(int gas_constituent_id)
{
	// binding to atmostracers
	return spec_heat_capacities_p_gas_lookup(gas_constituent_id);
}

double specific_gas_constants(int gas_constituent_id)
{
	// binding to atmostracers
	return specific_gas_constants_lookup(gas_constituent_id);
}






