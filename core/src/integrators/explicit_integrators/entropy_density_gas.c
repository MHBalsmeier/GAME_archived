/*
This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "../../spatial_operators/spatial_operators.h"

int integrate_entropy_density_gas(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{
	scalar_times_vector(state_old -> entropy_density_gas, state_new -> velocity_gas, diagnostics -> entropy_gas_flux_density, grid);
	divv_h(diagnostics -> entropy_gas_flux_density, forcings -> entropy_gas_flux_density_divv, grid);
	double rho_h, total_density;
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
	    if (config_info -> scalar_diffusion_on == 1)
	    {
	        total_density = state_old -> density_dry[i];
	        for (int k = 0; k < NO_OF_TRACERS; ++k)
	            total_density += state_old -> tracer_densities[k*NO_OF_SCALARS + i];
	        rho_h = state_old -> density_dry[i] + state_old -> tracer_densities[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i];
	        state_tendency -> entropy_density_gas[i] = -forcings -> entropy_gas_flux_density_divv[i] + 1/state_old -> temp_gas[i]*(rho_h/total_density*(diffusion_info -> temp_diffusion_heating[i] + config_info -> momentum_diffusion_on*diffusion_info -> heating_diss[i] + radiation_tendency[i]) + config_info -> tracers_on*config_info -> phase_transitions_on*diffusion_info -> tracer_heat_source_rates[NO_OF_CONDENSATED_TRACERS*NO_OF_SCALARS + i]);
	    }
	    else
	    {
	        state_tendency -> entropy_density_gas[i] = -forcings -> entropy_gas_flux_density_divv[i] + 1/state_old -> temp_gas[i]*radiation_tendency[i];
	    }
	}
	return 0;
}
