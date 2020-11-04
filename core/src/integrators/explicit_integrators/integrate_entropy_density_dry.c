/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "../../spatial_operators/spatial_operators.h"

int integrate_entropy_density_dry(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{
	// The dry specific entropy comes from the new time step.
	if (no_rk_step > 0)
	{
		for(int i = 0; i < NO_OF_SCALARS; ++i)
		{
			diagnostics -> specific_entropy_dry[i] = state_new -> entropy_density_dry[i]/state_new -> density_dry[i];
		}
	}
	// The dry mass flux density has already been updated.
	scalar_times_vector(diagnostics -> specific_entropy_dry, diagnostics -> mass_dry_flux_density, diagnostics -> entropy_dry_flux_density, grid, 0);
	divv_h(diagnostics -> entropy_dry_flux_density, forcings -> entropy_dry_flux_density_divv, grid);
	double rho_h, total_density;
	int k;
	#pragma omp parallel for private(rho_h, total_density, k)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
        total_density = state_old -> density_dry[i];
        for (k = 0; k < NO_OF_TRACERS; ++k)
        {
            total_density += state_old -> tracer_densities[k*NO_OF_SCALARS + i];
        }
        rho_h = state_old -> density_dry[i] + state_old -> tracer_densities[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i];
        state_tendency -> entropy_density_dry[i] =
        // convergence of adiabatic flux density
        -forcings -> entropy_dry_flux_density_divv[i]
        // diabatic source terms, need to be divided by the temperature, source terms have SI unit W/m^3
        // The source terms which act on the whole air must get the prefactor (dry density)/(total air density). 
         + 1/state_old -> temperature_gas[i]*(state_old -> density_dry[i]/total_density*(
         // This is the source term from the temperature diffusion.
         diffusion_info -> temperature_diffusion_heating[i]
         // This is the source term from the dissipation rate.
         + diffusion_info -> heating_diss[i]
         // This is the source term from radiative forcing.
         + radiation_tendency[i])
         // This is the power density arising from the phase transitions. They act on the gas phase and thus must get the prefactor (dry density)/(moist air density).
         + state_old -> density_dry[i]/rho_h*config_info -> phase_transitions_on*diffusion_info -> tracer_heat_source_rates[NO_OF_CONDENSED_TRACERS*NO_OF_SCALARS + i]);
	}
	return 0;
}
