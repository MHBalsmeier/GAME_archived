/*
This source file is part of the General Geophysical Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/game
*/

#include "../../enum_and_typedefs.h"
#include "../../spatial_operators/spatial_operators.h"
#include "../../diagnostics/diagnostics.h"

int integrate_continuity_dry(State *state_old, State *state_new, Interpolate_info *interpolation, State *state_tendency, Grid *grid, Dualgrid *dualgrid, Scalar_field radiation_tendency, Diagnostics *diagnostics, Forcings *forcings, Diffusion_info *diffusion_info, Config_info *config_info, int no_rk_step)
{
	// The dry mass flux density needs to be updated, becasue now, the new velocity value is known. This is a result of the forward-backward scheme.
	if (no_rk_step == 0)
	{
    	scalar_times_vector(state_old -> density_dry, state_new -> velocity_gas, diagnostics -> mass_dry_flux_density, grid);
	}
	else
	{
    	scalar_times_vector(state_new -> density_dry, state_new -> velocity_gas, diagnostics -> mass_dry_flux_density, grid);
	}
    divv_h(diagnostics -> mass_dry_flux_density, forcings -> mass_dry_flux_density_divv, grid);
    // This is where mass diffusion is updated.
    if (no_rk_step == 2)
    {
        grad(state_old -> density_dry, diffusion_info -> mass_dry_diffusion_flux_density, grid);
        calc_mass_diffusion_coeffs(config_info, state_old -> temp_gas, state_old -> density_dry, diffusion_info -> diffusion_coeff_numerical_h, diffusion_info -> diffusion_coeff_numerical_v);
        scalar_times_vector_scalar_h_v(diffusion_info -> diffusion_coeff_numerical_h, diffusion_info -> diffusion_coeff_numerical_v, diffusion_info -> mass_dry_diffusion_flux_density, diffusion_info -> mass_dry_diffusion_flux_density, grid);
        divv_h(diffusion_info -> mass_dry_diffusion_flux_density, diffusion_info -> mass_dry_diffusion_source_rate, grid);
    }
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
	    state_tendency -> density_dry[i] = -forcings -> mass_dry_flux_density_divv[i] + diffusion_info -> mass_dry_diffusion_source_rate[i];
    }
	return 0;
}
