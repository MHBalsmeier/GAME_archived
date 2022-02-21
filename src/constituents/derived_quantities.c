/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains functions calculating derived thermodynamic quantities of the atmosphere.
*/

#include <math.h>
#include "../game_constants.h"
#include "../game_types.h"
#include "constituents.h"

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
	/*
	This function calculates the specific heat capacity of the air at constant volume.
	*/
	
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
	/*
	This function calculates the specific heat capacity of the air at constant pressure.
	*/
	
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
	/*
	This function calculates the specific gas constant of the gas phase.
	*/
	
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
	/*
	This function calculates the density of the air.
	*/
	
	double result = 0;
	for (int i = 0; i < NO_OF_CONSTITUENTS; ++i)
	{
		result += state -> rho[i*NO_OF_SCALARS + grid_point_index];
	}
	return result;
}

double density_gas(State *state, int grid_point_index)
{
	/*
	This function calculates the density of the gas phase.
	*/
	
	double result = 0;
	for (int i = 0; i < NO_OF_GASEOUS_CONSTITUENTS; ++i)
	{
		result += state -> rho[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index];
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
    double thermal_velocity = sqrt(8.0*K_B*temperature/(M_PI*particle_mass));
    double particle_density = density/particle_mass;
    double cross_section = 4.0*M_PI*pow(particle_radius, 2.0);
    double mean_free_path = 1.0/(sqrt(2.0)*particle_density*cross_section);
    double result = 1.0/3.0*thermal_velocity*mean_free_path;
    return result;
}








